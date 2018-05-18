#include <assert.h>
#include "hkpriv.h"

/*********
 * Pairs *
 *********/

#include "ksort.h"
#define pair_lt(a, b) ((a).chr < (b).chr || ((a).chr == (b).chr && (a).pos < (b).pos))
KSORT_INIT(pair, struct hk_pair, pair_lt)

void hk_pair_sort(int32_t n_pairs, struct hk_pair *pairs)
{
	ks_introsort_pair(n_pairs, pairs);
}

int hk_pair_is_sorted(int32_t n_pairs, const struct hk_pair *pairs)
{
	int32_t i;
	for (i = 1; i < n_pairs; ++i) {
		const struct hk_pair *p, *q;
		p = &pairs[i-1], q = &pairs[i];
		if (p->chr < q->chr || (p->chr == q->chr && p->pos <= q->pos))
			continue;
		break;
	}
	return (i == n_pairs);
}

struct hk_pair *hk_seg2pair(int32_t n_segs, const struct hk_seg *segs, int min_dist, int max_seg, int min_mapq, int32_t *n_pairs_)
{
	int32_t m_pairs = 0, n_pairs = 0, i, j, st;
	struct hk_pair *pairs = 0;
	for (st = 0, i = 1; i <= n_segs; ++i) {
		if (i == n_segs || segs[i].frag_id != segs[i-1].frag_id) {
			if (i - st <= max_seg) {
				for (j = st + 1; j < i; ++j) {
					const struct hk_seg *s = &segs[j], *t = &segs[j-1];
					struct hk_pair *p;
					if (s->mapq < min_mapq || t->mapq < min_mapq)
						continue; // mapping quality too low
					if (n_pairs == m_pairs)
						EXPAND(pairs, m_pairs);
					p = &pairs[n_pairs++];
					if (t->chr < s->chr || (t->chr == s->chr && t->en < s->st)) {
						p->chr = (uint64_t)t->chr << 32 | s->chr;
						p->pos = (uint64_t)t->en  << 32 | s->st;
						p->phase[0]  = t->phase,  p->phase[1]  = s->phase;
						p->strand[0] = t->strand, p->strand[1] = s->strand;
					} else {
						p->chr = (uint64_t)s->chr << 32 | t->chr;
						p->pos = (uint64_t)s->st  << 32 | t->en;
						p->phase[0]  = s->phase,  p->phase[1]  = t->phase;
						p->strand[0] = s->strand, p->strand[1] = t->strand;
					}
					p->n = 0, p->tad_masked = 0;
					if (p->strand[0] * p->strand[1] >= 0 && t->chr == s->chr && (int32_t)p->pos - (int32_t)(p->pos>>32) < min_dist) {
						--n_pairs;
						continue;
					}
				}
			}
			st = i;
		}
	}
	ks_introsort_pair(n_pairs, pairs);
	*n_pairs_ = n_pairs;
	return pairs;
}

int32_t hk_pair_dedup(int n_pairs, struct hk_pair *pairs, int min_dist)
{
	int32_t i, j, n;
	for (i = n = 1; i < n_pairs; ++i) {
		struct hk_pair *q = &pairs[i];
		int32_t to_skip = 0;
		for (j = i - 1; j >= 0; --j) {
			struct hk_pair *p = &pairs[j];
			int32_t d;
			if (p->chr != q->chr || hk_ppos1(q) - hk_ppos1(p) >= min_dist)
				break;
			d = hk_ppos2(q) - hk_ppos2(p);
			if (d < min_dist && d > -min_dist && q->strand[0] == p->strand[0] && q->strand[1] == p->strand[1]) {
				to_skip = 1;
				break;
			}
		}
		if (!to_skip) pairs[n++] = pairs[i];
	}
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] duplicate rate: %.2f%% = %d / %d\n", __func__,
				100.0 * (n_pairs - n) / n_pairs, n_pairs - n, n_pairs);
	return n;
}

int32_t hk_pair_select_phased(int n_pairs, struct hk_pair *pairs)
{
	int32_t i, n;
	for (i = n = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (p->phase[0] < 0 && p->phase[1] < 0) continue;
		pairs[n++] = pairs[i];
	}
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] %d pairs remain\n", __func__, n);
	return n;
}

void hk_mask_by_tad(int32_t n_tads, const struct hk_pair *tads, int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t i, j, n_a, n_masked = 0;
	struct hk_pair *a;
	// populate _a_
	a = MALLOC(struct hk_pair, n_tads * 3);
	for (i = 0i, n_a = 0; i < n_tads; ++i) {
		const struct hk_pair *t = &tads[i];
		struct hk_pair *p;
		int32_t t1 = hk_ppos1(t), t2 = hk_ppos2(t);
		int32_t last_t2 = i > 0 && t->chr == (t-1)->chr? hk_ppos2(t-1) : 0;
		p = &a[n_a++], p->chr = t->chr, p->pos = (uint64_t)last_t2 << 32 | t1;
		p = &a[n_a++], p->chr = t->chr, p->pos = (uint64_t)t1 << 32 | t2;
		if (i == n_tads - 1 || t->chr != (t+1)->chr)
			p = &a[n_a++], p->chr = t->chr, p->pos = (uint64_t)t2 << 32 | INT32_MAX;
	}
	// filter
	for (i = j = n_masked = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		int kept = 1;
		if (p->chr>>32 == (int32_t)p->chr) { // intra-chromosomal pairs
			const struct hk_pair *t;
			int32_t p1 = hk_ppos1(p);
			while (j < n_a && (a[j].chr < p->chr || (a[j].chr == p->chr && hk_ppos2(&a[j]) <= p1)))
				++j;
			if (j == n_a) break;
			t = &a[j];
			if (p->chr == t->chr && hk_ppos1(t) <= p1 && hk_ppos2(p) <= hk_ppos2(t)) // contained in TAD
				kept = 0;
		}
		pairs[i].tad_masked = !kept;
		if (!kept) ++n_masked;
	}
	free(a);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] masked %d out of %d pairs\n", __func__, n_masked, n_pairs);
}

struct hk_map *hk_pair_sep_phase(const struct hk_map *m, float phase_thres)
{
	int32_t i, m_ppairs = 0, *ploidy_XY, *old2new;
	struct hk_map *p;

	ploidy_XY = hk_sd_ploidy_XY(m->d, 0);
	old2new = CALLOC(int32_t, m->d->n);
	for (i = 1; i < m->d->n; ++i)
		old2new[i] = old2new[i-1] + (ploidy_XY[i-1]>>8);

	p = CALLOC(struct hk_map, 1);
	p->d = hk_sd_sep_phase(m->d, ploidy_XY);
	for (i = 0; i < m->n_pairs; ++i) {
		const struct hk_pair *q = &m->pairs[i];
		struct hk_pair *r;
		int32_t j, max_j = -1, chr[2];
		float max = -1e30f;
		for (j = 0; j < 4; ++j)
			if (q->_.p4[j] > max)
				max = q->_.p4[j], max_j = j;
		if (max < phase_thres) continue;
		chr[0] = q->chr>>32, chr[1] = (int32_t)q->chr;
		chr[0] = old2new[chr[0]] + (ploidy_XY[chr[0]]>>8 == 1? 0 : max_j>>1&1);
		chr[1] = old2new[chr[1]] + (ploidy_XY[chr[1]]>>8 == 1? 0 : max_j&1);
		if (p->n_pairs == m_ppairs)
			EXPAND(p->pairs, m_ppairs);
		r = &p->pairs[p->n_pairs++];
		*r = *q;
		if (chr[0] > chr[1]) {
			r->chr = (uint64_t)chr[1] << 32 | chr[0];
			r->pos = q->pos<<32 | q->pos>>32;
			r->phase[0]  = q->phase[1],  r->phase[1]  = q->phase[0];
			r->strand[0] = q->strand[1], r->strand[1] = q->strand[0];
		} else r->chr = (uint64_t)chr[0] << 32 | chr[1];
		r->_.phased_prob = max;
	}
	free(old2new);
	free(ploidy_XY);
	hk_pair_sort(p->n_pairs, p->pairs);
	return p;
}

int32_t hk_pair_filter(int32_t n_pairs, struct hk_pair *pairs, int32_t max_radius, int32_t min_cnt, float drop_frac)
{
	int32_t i, k, *cnt, min;
	hk_pair_count_nei(n_pairs, pairs, max_radius);
	cnt = CALLOC(int32_t, n_pairs);
	for (i = 0; i < n_pairs; ++i)
		cnt[i] = pairs[i].n;
	min = ks_ksmall_int32_t(n_pairs, cnt, n_pairs * drop_frac);
	min = min > min_cnt? min : min_cnt;
	free(cnt);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] threshold: %d\n", __func__, min);
	for (i = k = 0; i < n_pairs; ++i)
		if (pairs[i].n >= min)
			pairs[k++] = pairs[i];
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] filtered out %d isolated pairs\n", __func__, n_pairs - k);
	return k;
}

/********************
 * Print segs/pairs *
 ********************/

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
	if (flag & HK_OUT_CNT) fprintf(fp, " count");
	if (flag & HK_OUT_P4) fprintf(fp, " p00 p01 p10 p11");
	if (flag & HK_OUT_PHASE) fprintf(fp, " phase1 phase2");
	fputc('\n', fp);
	for (i = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		fprintf(fp, ".\t%s\t%d\t%s\t%d\t%c\t%c", d->name[p->chr>>32], (int32_t)(p->pos>>32),
				d->name[(int32_t)p->chr], (int32_t)p->pos,
				p->strand[0] > 0? '+' : p->strand[0] < 0? '-' : '.',
				p->strand[1] > 0? '+' : p->strand[1] < 0? '-' : '.');
		if (flag & HK_OUT_CNT) fprintf(fp, "\t%d", p->n);
		if (flag & HK_OUT_P4) fprintf(fp, "\t%.3f\t%.3f\t%.3f\t%.3f", p->_.p4[0], p->_.p4[1], p->_.p4[2], p->_.p4[3]);
		if (flag & HK_OUT_PHASE) fprintf(fp, "\t%c\t%c", p->phase[0] < 0? '.' : '0' + p->phase[0], p->phase[1] < 0? '.' : '0' + p->phase[1]);
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
		fprintf(fp, "%s\t%d\t%f\t%f\t%f\n", m->d->name[b->chr], b->st,
				m->x[i][0], m->x[i][1], m->x[i][2]);
	}
}
