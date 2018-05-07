#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "hickit.h"
#include "hkpriv.h"

/**********************
 * Default parameters *
 **********************/

int hk_verbose = 3;

void hk_opt_init(struct hk_opt *c)
{
	memset(c, 0, sizeof(struct hk_opt));
	c->min_dist = 1000;
	c->max_seg = 3;
	c->min_mapq = 20;
	c->min_flt_cnt = 0;
	c->min_tad_size = 10;
	c->area_weight = 5.0f;
	c->min_radius = 50000;
	c->max_radius = 10000000;
	c->max_nei = 50;
	c->pseudo_cnt = 0.4f;
	c->n_iter = 1000;
	c->n_burnin = 1000;
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

struct hk_sdict *hk_sd_sep_phase(const struct hk_sdict *d, int32_t *ploidy_XY)
{
	struct hk_sdict *dp;
	int i;
	dp = hk_sd_init();
	for (i = 0; i < d->n; ++i) {
		int ploidy = ploidy_XY? ploidy_XY[i] >> 8 : 1;
		if (ploidy <= 1) {
			hk_sd_put(dp, d->name[i], d->len[i]);
		} else {
			int l, j;
			char *s;
			l = strlen(d->name[i]);
			s = CALLOC(char, l + 2);
			strcpy(s, d->name[i]);
			s[l+1] = 0;
			for (j = 0; j < ploidy; ++j) {
				s[l] = 'a' + j;
				hk_sd_put(dp, s, d->len[i]);
			}
			free(s);
		}
	}
	return dp;
}

struct hk_sdict *hk_sd_dup(const struct hk_sdict *d)
{
	return hk_sd_sep_phase(d, 0);
}

int32_t *hk_sd_ploidy_XY(const struct hk_sdict *d, int32_t *sex_flag)
{
	int32_t i, *ploidy, n_X = 0, n_Y = 0;
	if (sex_flag) *sex_flag = 0;
	for (i = 0; i < d->n; ++i) {
		int32_t is_X, is_Y;
		is_X = (strchr(d->name[i], 'X') != 0);
		is_Y = (strchr(d->name[i], 'Y') != 0);
		if (is_X && is_Y) continue;
		if (is_X) ++n_X;
		if (is_Y) ++n_Y;
	}
	if (n_X == 0 && n_Y == 0)
		fprintf(stderr, "[W::%s] no sex chromosomes (identified by 'X' or 'Y' in chr names)\n", __func__);
	if (n_X > 1 || n_Y > 1) {
		fprintf(stderr, "[E::%s] multiple chr contain 'X' or 'Y' in names\n", __func__);
		exit(1);
	}
	if (n_X > 0 && sex_flag) *sex_flag |= 1;
	if (n_Y > 0 && sex_flag) *sex_flag |= 2;
	ploidy = CALLOC(int32_t, d->n);
	for (i = 0; i < d->n; ++i) {
		int32_t is_X, is_Y;
		is_X = (strchr(d->name[i], 'X') != 0);
		is_Y = (strchr(d->name[i], 'Y') != 0);
		if (is_X && is_Y) continue;
		if (n_Y == 0) ploidy[i] = 2 << 8;
		else ploidy[i] = (is_X || is_Y? 1 : 2) << 8;
		if (is_X) ploidy[i] |= 1;
		if (is_Y) ploidy[i] |= 2;
	}
	return ploidy;
}

/*************************
 * Contact alloc/dealloc *
 *************************/

struct hk_map *hk_map_init(void)
{
	struct hk_map *m;
	m = CALLOC(struct hk_map, 1);
	m->d = hk_sd_init();
	return m;
}

void hk_map_destroy(struct hk_map *m)
{
	free(m->pairs);
	free(m->segs);
	hk_sd_destroy(m->d);
	free(m);
}

/*************************
 * Seg/pairs file parser *
 *************************/

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

static int64_t hk_parse_64(const char *s, char **t, int *has_digit)
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

static void hk_parse_seg(struct hk_seg *s, struct hk_sdict *d, int32_t frag_id, char *str)
{
	char *p, *q;
	int i, has_digit;

	s->frag_id = frag_id;
	s->chr = s->st = s->en = -1;
	s->strand = 0, s->phase = -1, s->mapq = 0;

	for (i = 0, p = q = str;; ++q) {
		if (*q == HK_SUB_DELIM || *q == 0) {
			int c = *q;
			if (i == 0) {
				assert(c != 0);
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
			if (c == 0) break;
			p = q + 1;
		}
	}
}

static void hk_parse_pair(struct hk_pair *p, struct hk_sdict *d, int n_fields, char **fields)
{
	int32_t c1, c2;
	int64_t p1, p2;
	char *q;
	int j, has_digit;
	c1 = hk_sd_put(d, fields[1], 0);
	p1 = hk_parse_64(fields[2], &q, &has_digit);
	assert(p1 >= 0 && has_digit);
	c2 = hk_sd_put(d, fields[3], 0);
	p2 = hk_parse_64(fields[4], &q, &has_digit);
	assert(p2 >= 0 && has_digit);
	memset(p, 0, sizeof(struct hk_pair));
	p->chr = (uint64_t)c1 << 32 | c2;
	p->pos = (uint64_t)p1 << 32 | p2;
	p->phase[0] = p->phase[1] = -1;
	if (n_fields >= 7) {
		p->strand[0] = *fields[5] == '+'? 1 : *fields[5] == '-'? -1 : 0;
		p->strand[1] = *fields[6] == '+'? 1 : *fields[6] == '-'? -1 : 0;
		if (n_fields >= 11) { // FIXME: make this more general
			for (j = 0; j < 4; ++j)
				p->_.p4[j] = atof(fields[7 + j]);
		} else if (n_fields >= 9) {
			p->phase[0] = *fields[7] == '.'? -1 : (int)*fields[7] - '0';
			p->phase[1] = *fields[8] == '.'? -1 : (int)*fields[8] - '0';
		}
	}
}

static inline int is_pos_int(const char *s)
{
	for (; *s; ++s)
		if (*s < '0' || *s > '9') return 0;
	return 1;
}

static void parse_chr(struct hk_sdict *d, char *s)
{
	char *chr, *p, *q;
	int64_t len;
	int has_digit;
	for (p = s + 12; isspace(*p) && *p != 0; ++p) {}
	assert(*p);
	for (q = p; *q != 0 && !isspace(*q); ++q) {}
	assert(*q);
	*q = 0, chr = p, p = q + 1;
	len = hk_parse_64(p, &p, &has_digit);
	assert(has_digit && len > 0 && len <= INT32_MAX);
	hk_sd_put(d, chr, len);
}

static char **split_fields(char *s, char **fields, int32_t *n_fields_, int32_t *m_fields_)
{
	int32_t n_fields, m_fields = *m_fields_;
	char *p, *q;
	for (n_fields = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == 0) {
			int c = *q;
			if (n_fields == m_fields)
				EXPAND(fields, m_fields);
			fields[n_fields++] = p;
			*q = 0, p = q + 1;
			if (c == 0) break;
		}
	}
	*n_fields_ = n_fields, *m_fields_ = m_fields;
	return fields;
}

struct hk_map *hk_map_read(const char *fn)
{
	gzFile fp;
	kstring_t str = {0,0,0};
	kstream_t *ks;
	int dret;
	int32_t m_segs = 0, m_pairs = 0, n_fields = 0, m_fields = 0;
	int64_t n_data_rows = 0;
	char **fields = 0;
	struct hk_map *m;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	if (fp == 0) return 0;
	m = hk_map_init();
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		int32_t k, n_segs = 0;
		if (str.l && str.s[0] != '#') ++n_data_rows;
		if (n_data_rows == 0 && str.l >= 12 + 3 && strncmp(str.s, "#chromosome:", 12) == 0)
			parse_chr(m->d, str.s);
		if (str.s[0] == '#') continue;
		fields = split_fields(str.s, fields, &n_fields, &m_fields);
		if (n_fields < 5) goto parse_seg;
		if (!is_pos_int(fields[2]) || !is_pos_int(fields[4])) goto parse_seg;
		if (n_fields >= 7 && (strlen(fields[5]) != 1 || strlen(fields[6]) != 1))
			goto parse_seg;
		if (m->n_pairs == m_pairs)
			EXPAND(m->pairs, m_pairs);
		hk_parse_pair(&m->pairs[m->n_pairs++], m->d, n_fields, fields);
		continue;
parse_seg:
		for (k = 1; k < n_fields; ++k) {
			if (m->n_segs == m_segs)
				EXPAND(m->segs, m_segs);
			hk_parse_seg(&m->segs[m->n_segs++], m->d, m->n_frags, fields[k]);
			++n_segs;
		}
		if (n_segs > 0) ++m->n_frags;
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	if (m->pairs && !hk_pair_is_sorted(m->n_pairs, m->pairs))
		hk_pair_sort(m->n_pairs, m->pairs);
	return m;
}

struct hk_bmap *hk_3dg_read(const char *fn)
{
	gzFile fp;
	kstring_t str = {0,0,0};
	kstream_t *ks;
	int32_t i, dret, m_beads = 0, n_fields = 0, m_fields = 0, n_data_rows = 0;
	char **fields = 0;
	struct hk_bmap *m;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	if (fp == 0) return 0;
	m = CALLOC(struct hk_bmap, 1);
	m->d = hk_sd_init();
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		if (str.l && str.s[0] != '#') ++n_data_rows;
		if (n_data_rows == 0 && str.l >= 12 + 3 && strncmp(str.s, "#chromosome:", 12) == 0)
			parse_chr(m->d, str.s);
		if (str.s[0] == '#') continue;
		fields = split_fields(str.s, fields, &n_fields, &m_fields);
		if (n_fields >= 5) {
			struct hk_bead *b;
			if (m->n_beads == m_beads) {
				EXPAND(m->beads, m_beads);
				REALLOC(m->x, m_beads);
			}
			m->x[m->n_beads][0] = atof(fields[2]);
			m->x[m->n_beads][1] = atof(fields[3]);
			m->x[m->n_beads][2] = atof(fields[4]);
			b = &m->beads[m->n_beads++];
			b->chr = hk_sd_put(m->d, fields[0], 0);
			b->st = b->en = atoi(fields[1]);
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	for (i = 1; i <= m->n_beads; ++i) {
		if (i == m->n_beads || m->beads[i].chr != m->beads[i-1].chr)
			m->beads[i-1].en = m->d->len[m->beads[i-1].chr];
		else m->beads[i-1].en = m->beads[i].st;
	}
	hk_bmap_set_offcnt(m);
	return m;
}

void hk_map_phase_male_XY(struct hk_map *m)
{
	int32_t i, sex_flag, *ploidy_XY;
	ploidy_XY = hk_sd_ploidy_XY(m->d, &sex_flag);
	if (sex_flag & 2) { // chrY present, a male
		if (m->pairs) {
			for (i = 0; i < m->n_pairs; ++i) {
				struct hk_pair *p = &m->pairs[i];
				int32_t chr[2];
				chr[0] = p->chr >> 32;
				chr[1] = (int32_t)p->chr;
				if (ploidy_XY[chr[0]]&1) p->phase[0] = 1;
				else if (ploidy_XY[chr[0]]&2) p->phase[0] = 0;
				if (ploidy_XY[chr[1]]&1) p->phase[1] = 1;
				else if (ploidy_XY[chr[1]]&2) p->phase[1] = 0;
			}
		}
		if (m->segs) {
			for (i = 0; i < m->n_segs; ++i) {
				struct hk_seg *p = &m->segs[i];
				if (ploidy_XY[p->chr]&1) p->phase = 1;
				else if (ploidy_XY[p->chr]&2) p->phase = 0;
			}
		}
	}
	free(ploidy_XY);
}
