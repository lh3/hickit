#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "hickit.h"
#include "hkpriv.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

/*************
 * Utilities *
 *************/

int64_t hk_parse_64(const char *s, char **t, int *has_digit)
{
	const char *p = s;
	int is_neg = 0;
	int64_t x = 0;
	*has_digit = 0;
	while (isspace(*p)) ++p;
	if (*p == '-') is_neg = 1, ++p;
	else if (*p == '+') ++p;
	while (*p >= '0' && *p <= '9') {
		x = x * 10 + (*p - '0');
		*has_digit = 1;
		++p;
	}
	if (is_neg) x = -x;
	*t = (char*)p;
	return x;
}

/*********************
 * String dictionary *
 *********************/

#include "khash.h"
KHASH_MAP_INIT_STR(str, int32_t)
typedef khash_t(str) sdict_t;

struct hk_sdict *hk_sd_init(void)
{
	struct hk_sdict *d;
	d = CALLOC(struct hk_sdict, 1);
	d->h = kh_init(str);
	return d;
}

void hk_sd_destroy(struct hk_sdict *d)
{
	int32_t i;
	for (i = 0; i < d->n; ++i)
		free(d->name[i]);
	free(d->name);
	free(d->len);
	kh_destroy(str, (sdict_t*)d->h);
}

int32_t hk_sd_put(struct hk_sdict *d, const char *s, int32_t len)
{
	sdict_t *h = (sdict_t*)d->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, s, &absent);
	if (absent) {
		if (d->n == d->m) {
			EXPAND(d->name, d->m);
			REALLOC(d->len, d->m);
		}
		d->len[d->n] = len;
		kh_key(h, k) = d->name[d->n] = strdup(s);
		kh_val(h, k) = d->n++;
	}
	return kh_val(h, k);
}

/***************
 * Segment I/O *
 ***************/

struct hk_map *hk_map_init(void)
{
	struct hk_map *m;
	m = CALLOC(struct hk_map, 1);
	m->d = hk_sd_init();
	return m;
}

void hk_map_destroy(struct hk_map *m)
{
	hk_sd_destroy(m->d);
	free(m->seg);
	free(m);
}

void hk_parse_seg(struct hk_seg *s, struct hk_sdict *d, int32_t frag_id, char *start, char *end)
{
	char *p, *q;
	int i, has_digit;

	s->frag_id = frag_id;
	s->chr = s->st = s->en = -1;
	s->strand = 0, s->phase = -1, s->mapq = 0;

	for (i = 0, p = q = start;; ++q) {
		if (q == end || *q == HK_SUB_DELIM) {
			if (i == 0) {
				assert(q < end);
				*q = 0;
				s->chr = hk_sd_put(d, p, 0);
				*q = HK_SUB_DELIM;
			} else if (i == 1) {
				s->st = s->en = hk_parse_64(p, &p, &has_digit);
				assert(has_digit && s->st >= 0);
			} else if (i == 2) {
				s->en = hk_parse_64(p, &p, &has_digit);
				if (!has_digit || s->en < s->st) s->en = s->st;
			} else if (i == 3) {
				if (*p == '+') s->strand = 1;
				else if (*p == '-') s->strand = -1;
			} else if (i == 4) {
				if (*p == '0') s->phase = 0;
				else if (*p == '1') s->phase = 1;
			} else if (i == 5) {
				s->mapq = hk_parse_64(p, &p, &has_digit);
				if (!has_digit || s->mapq < 0 || s->mapq >= 255) s->mapq = 0;
			}
			++i;
			if (q == end) break;
			p = q + 1;
		}
	}
}

struct hk_map *hk_map_read(const char *fn)
{
	gzFile fp;
	kstring_t str = {0,0,0};
	kstream_t *ks;
	int dret;
	struct hk_map *m;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	if (fp == 0) return 0;
	m = hk_map_init();
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		char *p, *q;
		int32_t k, n_seg = 0;
		if (str.l >= 12 + 3 && strncmp(str.s, "#chromosome:", 12) == 0) {
			char *chr;
			int64_t len;
			int has_digit;
			for (p = str.s + 12; isspace(*p) && *p != 0; ++p) {}
			assert(*p);
			for (q = p; *q != 0 && !isspace(*q); ++q) {}
			assert(*q);
			*q = 0, chr = p, p = q + 1;
			len = hk_parse_64(p, &p, &has_digit);
			assert(has_digit && len > 0 && len <= INT32_MAX);
			hk_sd_put(m->d, chr, len);
		}
		for (k = 0, p = q = str.s;; ++q) {
			if (*q == '\t' || *q == 0) {
				if (k > 0) {
					if (m->n_seg == m->m_seg)
						EXPAND(m->seg, m->m_seg);
					hk_parse_seg(&m->seg[m->n_seg++], m->d, m->n_frag, p, q);
					++n_seg;
				}
				++k, p = q + 1;
				if (*q == 0) break;
			}
		}
		if (n_seg > 0) ++m->n_frag;
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return m;
}

void hk_print_chr(FILE *fp, const struct hk_map *m)
{
	int32_t i;
	for (i = 0; i < m->d->n; ++i)
		if (m->d->len[i] > 0)
			fprintf(fp, "#chromosome: %s %d\n", m->d->name[i], m->d->len[i]);
}

void hk_map_print(FILE *fp, const struct hk_map *m)
{
	int32_t i, last_frag = -1;
	hk_print_chr(fp, m);
	for (i = 0; i < m->n_seg; ++i) {
		struct hk_seg *s = &m->seg[i];
		if (s->frag_id != last_frag) {
			if (last_frag >= 0) fputc('\n', fp);
			fputc('.', fp);
			last_frag = s->frag_id;
		}
		fprintf(fp, "\t%s%c%d%c", m->d->name[s->chr], HK_SUB_DELIM, s->st, HK_SUB_DELIM);
		if (s->en > s->st) fprintf(fp, "%d", s->en);
		fprintf(fp, "%c%c%c%c%c%d", HK_SUB_DELIM, s->strand > 0? '+' : s->strand < 0? '-' : '.',
				HK_SUB_DELIM, s->phase == 0? '0' : s->phase == 1? '1' : '.',
				HK_SUB_DELIM, s->mapq);
	}
	fputc('\n', fp);
}
