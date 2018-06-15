#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "hickit.h"
#include "hkpriv.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

char *hk_pair_cols[] = { // when modify this array, append; DON'T insert in the middle
	"phase0",        // 0
	"phase1",        // 1
	"phase_prob00",  // 2
	"phase_prob01",  // 3
	"phase_prob10",  // 4
	"phase_prob11",  // 5
	"n_neighbors",   // 6
	"n_contained",   // 7
	"prob",          // 8
	"n_nei_corner",  // 9
	"peak_density0", // 10
	"peak_density1", // 11
	"peak_density2", // 12
	NULL
};

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

/******************
 * Format parsers *
 ******************/

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

static void hk_parse_pair(struct hk_pair *p, struct hk_sdict *d, int n_extra_cols, int *extra_cols, int n_fields, char **fields)
{
	int32_t c1, c2;
	int64_t p1, p2;
	char *q;
	int i, has_digit;

	assert(n_fields == 5 || n_fields >= 7); // TODO: turn this into human readable errors
	if (n_fields >= 7)
		assert(n_extra_cols == n_fields - 7);

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
		for (i = 0; i < n_extra_cols; ++i) {
			int c = i + 7, e = extra_cols[i];
			if (e == 0 || e == 1) { // phase0 or phase1
				p->phase[e] = *fields[c] == '.'? -1 : (int)*fields[c] - '0';
			} else if (e >= 2 && e <= 5) { // phase_prob00, 01, 10 and 11
				p->_.p4[e - 2] = atof(fields[c]);
			} else if (e == 6) { // n_neighbors
				p->n_nei = atoi(fields[c]);
			} else if (e == 7) { // phased_prob
				p->_.phased_prob = atof(fields[c]);
			} else if (e == 8) { // n_contained
				p->n_ctn = atoi(fields[c]);
			} else if (e == 9) { // n_nei_corner
				p->n_nei_corner = atoi(fields[c]);
			} else if (e >= 10 && e <= 12) { // peak_density0, 1 and 2
				p->_.peak_density[e - 10] = atof(fields[c]);
			}
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

static int *parse_pair_cols(char *s, int *n_extra_cols_)
{
	char *p, *q;
	int i, k, n_cols, n_extra_cols, *extra_cols = 0;

	for (p = s + 9; isspace(*p) && *p != 0; ++p) {} // skip spaces following ':'
	if (p == 0) return 0;
	for (q = p, n_cols = 0; *q; ++q)
		if (isspace(*q)) ++n_cols;
	++n_cols;
	if (n_cols <= 7) return 0; // no custom columns

	*n_extra_cols_ = n_extra_cols = n_cols - 7;
	extra_cols = CALLOC(int, n_extra_cols);
	for (q = p, k = 0;; ++q) {
		if (*q == 0 || isspace(*q)) {
			if (k >= 7) {
				for (i = 0; hk_pair_cols[i]; ++i)
					if (strncmp(p, hk_pair_cols[i], q - p) == 0)
						break;
				if (hk_pair_cols[i] == 0 && hk_verbose >= 2) {
					int c = *q;
					*q = 0;
					fprintf(stderr, "[W::%s] unrecognized column \"%s\" will be ignored\n", __func__, p);
					*q = c;
				}
				extra_cols[k - 7] = hk_pair_cols[i]? i : -1;
			}
			++k;
			if (*q == 0) break;
			p = q + 1;
		}
	}
	return extra_cols;
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

void hk_map_set_cols(struct hk_map *m, int n_extra_cols, int *extra_cols)
{
	m->cols = 0;
	if (m->segs) {
		int i, t[2];
		for (i = 0, t[0] = t[1] = 0; i < m->n_segs; ++i) {
			if (m->segs[i].phase >= 0 && m->segs[i].phase <= 1)
				t[m->segs[i].phase] = 1;
			if (t[0] && t[1]) break;
		}
		if (i < m->n_segs) m->cols |= 3;
	} else if (m->pairs) {
		int i, extra_flags, t[2];
		for (i = 0, extra_flags = 0; i < n_extra_cols; ++i)
			extra_flags |= 1 << extra_cols[i];
		if ((extra_flags & 0x3c) == 0x3c)
			m->cols |= 0x3c;
		for (i = 0, t[0] = t[1] = 0; i < m->n_pairs; ++i) {
			if (m->pairs[i].phase[0] >= 0) t[0] = 1;
			if (m->pairs[i].phase[1] >= 0) t[1] = 1;
			if (t[0] && t[1]) break;
		}
		if (i < m->n_pairs) m->cols |= 3;
	}
}

struct hk_map *hk_map_read(const char *fn)
{
	gzFile fp;
	kstring_t str = {0,0,0};
	kstream_t *ks;
	int dret, *extra_cols = 0, n_extra_cols = 0, is_pairs = 0;
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
		if (n_data_rows == 0) {
			if (str.l >= 12 + 3 && strncmp(str.s, "#chromosome:", 12) == 0)
				parse_chr(m->d, str.s);
			else if (str.l >= 9 && strncmp(str.s, "#columns:", 9) == 0) {
				extra_cols = parse_pair_cols(str.s, &n_extra_cols);
				is_pairs = 1;
			}
		}
		if (str.s[0] == '#') continue;
		fields = split_fields(str.s, fields, &n_fields, &m_fields);
		if (!is_pairs) {
			if (n_fields < 5) goto parse_seg;
			if (!is_pos_int(fields[2]) || !is_pos_int(fields[4])) goto parse_seg;
			if (n_fields >= 7 && (strlen(fields[5]) != 1 || strlen(fields[6]) != 1))
				goto parse_seg;
		}
		if (m->n_pairs == m_pairs)
			EXPAND(m->pairs, m_pairs);
		hk_parse_pair(&m->pairs[m->n_pairs++], m->d, n_extra_cols, extra_cols, n_fields, fields);
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
	hk_map_set_cols(m, n_extra_cols, extra_cols);
	free(fields);
	free(extra_cols);
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);

	if (m->pairs && !hk_pair_is_sorted(m->n_pairs, m->pairs))
		hk_pair_sort(m->n_pairs, m->pairs);
	if (hk_verbose >= 3) {
		if (m->n_segs) fprintf(stderr, "[M::%s] read %d segments\n", __func__, m->n_segs);
		if (m->n_pairs) fprintf(stderr, "[M::%s] read %d pairs\n", __func__, m->n_pairs);
	}
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
			if (n_fields >= 6) {
				REALLOC(m->feat, m_beads);
				m->feat[m->n_beads] = atof(fields[5]);
			}
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

/*****************
 * Format output *
 *****************/

static void hk_print_chr(FILE *fp, const struct hk_sdict *d)
{
	int32_t i;
	for (i = 0; i < d->n; ++i)
		if (d->len[i] > 0)
			fprintf(fp, "#chromosome: %s %d\n", d->name[i], d->len[i]);
}

void hk_print_seg(FILE *fp, const struct hk_sdict *d, int32_t n_segs, const struct hk_seg *segs)
{
	int32_t i, last_frag = -1;
	hk_print_chr(fp, d);
	for (i = 0; i < n_segs; ++i) {
		const struct hk_seg *s = &segs[i];
		if (s->frag_id != last_frag) {
			if (last_frag >= 0) fputc('\n', fp);
			fputc('.', fp);
			last_frag = s->frag_id;
		}
		fprintf(fp, "\t%s%c%d%c", d->name[s->chr], HK_SUB_DELIM, s->st, HK_SUB_DELIM);
		if (s->en > s->st) fprintf(fp, "%d", s->en);
		fprintf(fp, "%c%c%c%c%c%d", HK_SUB_DELIM, s->strand > 0? '+' : s->strand < 0? '-' : '.',
				HK_SUB_DELIM, s->phase == 0? '0' : s->phase == 1? '1' : '.',
				HK_SUB_DELIM, s->mapq);
	}
	fputc('\n', fp);
}

void hk_print_pair(FILE *fp, int flag, const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs)
{
	int32_t i;
	fprintf(fp, "## pairs format v1.0\n");
	fprintf(fp, "#sorted: chr1-chr2-pos1-pos2\n");
	fprintf(fp, "#shape: upper triangle\n");
	hk_print_chr(fp, d);
	fprintf(fp, "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2");
	if ((flag & 3) == 3) fprintf(fp, " %s %s", hk_pair_cols[0], hk_pair_cols[1]);
	if (flag & 1<<6) fprintf(fp, " %s", hk_pair_cols[6]);
	if (flag & 1<<9) fprintf(fp, " %s", hk_pair_cols[9]);
	if (flag & 1<<7) fprintf(fp, " %s", hk_pair_cols[7]);
	if ((flag & 0x3c) == 0x3c) fprintf(fp, " %s %s %s %s", hk_pair_cols[2], hk_pair_cols[3], hk_pair_cols[4], hk_pair_cols[5]);
	else if (flag & 1<<8) fprintf(fp, " %s", hk_pair_cols[8]);
	else if ((flag & 0x1c00) == 0x1c00) fprintf(fp, " %s %s %s", hk_pair_cols[10], hk_pair_cols[11], hk_pair_cols[12]);
	fputc('\n', fp);
	for (i = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		fprintf(fp, ".\t%s\t%d\t%s\t%d\t%c\t%c", d->name[p->chr>>32], (int32_t)(p->pos>>32),
				d->name[(int32_t)p->chr], (int32_t)p->pos,
				p->strand[0] > 0? '+' : p->strand[0] < 0? '-' : '.',
				p->strand[1] > 0? '+' : p->strand[1] < 0? '-' : '.');
		if ((flag & 3) == 3) fprintf(fp, "\t%c\t%c", p->phase[0] < 0? '.' : '0' + p->phase[0], p->phase[1] < 0? '.' : '0' + p->phase[1]);
		if (flag & 1<<6) fprintf(fp, "\t%d", p->n_nei);
		if (flag & 1<<9) fprintf(fp, "\t%d", p->n_nei_corner);
		if (flag & 1<<7) fprintf(fp, "\t%d", p->n_ctn);
		if ((flag & 0x3c) == 0x3c) fprintf(fp, "\t%.3f\t%.3f\t%.3f\t%.3f", p->_.p4[0], p->_.p4[1], p->_.p4[2], p->_.p4[3]);
		else if (flag & 1<<8) fprintf(fp, "\t%.4f", p->_.phased_prob);
		else if ((flag & 0x1c00) == 0x1c00) fprintf(fp, "\t%.4f\t%.4f\t%.4f", p->_.peak_density[0]*1e6, p->_.peak_density[1]*1e6, p->_.peak_density[2]*1e6);
		fputc('\n', fp);
	}
}

void hk_print_bmap(FILE *fp, const struct hk_bmap *m)
{
	int32_t i;
	fprintf(fp, "## pairs format v1.0\n");
	fprintf(fp, "#sorted: chr1-pos1-chr2-pos2\n");
	fprintf(fp, "#shape: upper triangle\n");
	hk_print_chr(fp, m->d);
	for (i = 0; i < m->n_pairs; ++i) {
		const struct hk_bpair *p = &m->pairs[i];
		const struct hk_bead *b[2];
		b[0] = &m->beads[p->bid[0]];
		b[1] = &m->beads[p->bid[1]];
		fprintf(fp, ".\t%s\t%d\t%s\t%d\t%d\t%.4f\n", m->d->name[b[0]->chr], b[0]->st,
				m->d->name[b[1]->chr], b[1]->st, p->n, p->p);
	}
}

void hk_print_3dg(FILE *fp, const struct hk_bmap *m)
{
	int32_t i;
	hk_print_chr(fp, m->d);
	for (i = 0; i < m->n_beads; ++i) {
		const struct hk_bead *b = &m->beads[i];
		fprintf(fp, "%s\t%d\t%f\t%f\t%f", m->d->name[b->chr], b->st,
				m->x[i][0], m->x[i][1], m->x[i][2]);
		if (m->feat) fprintf(fp, "\t%f", m->feat[i]);
		fputc('\n', fp);
	}
}
