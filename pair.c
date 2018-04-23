#include <assert.h>
#include "hkpriv.h"

/*********
 * Pairs *
 *********/

#include "ksort.h"
#define pair_lt(a, b) ((a).chr < (b).chr || ((a).chr == (b).chr && (a).pos < (b).pos))
KSORT_INIT(pair, struct hk_pair, pair_lt)

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
					p->n = 0;
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

int32_t hk_mask_by_tad(int32_t n_tads, const struct hk_pair *tads, int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t i, j, k, n_a;
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
	for (i = j = k = 0; i < n_pairs; ++i) {
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
		if (kept) pairs[k++] = pairs[i];
	}
	for (; i < n_pairs; ++i) // copy the rest over
		pairs[k++] = pairs[i];
	free(a);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] masked %d out of %d pairs\n", __func__, n_pairs - k, n_pairs);
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
	if (flag & HK_OUT_PHASE) fprintf(fp, " phase1 phase2");
	fputc('\n', fp);
	for (i = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		fprintf(fp, ".\t%s\t%d\t%s\t%d\t%c\t%c", d->name[p->chr>>32], (int32_t)(p->pos>>32),
				d->name[(int32_t)p->chr], (int32_t)p->pos,
				p->strand[0] > 0? '+' : p->strand[0] < 0? '-' : '.',
				p->strand[1] > 0? '+' : p->strand[1] < 0? '-' : '.');
		if (flag & HK_OUT_PHASE) {
			if (flag & HK_OUT_PHASE_REAL)
				fprintf(fp, "\t%.3f\t%.3f", p->_.phase_prob[0], p->_.phase_prob[1]);
			else
				fprintf(fp, "\t%c\t%c", p->phase[0] < 0? '.' : '0' + p->phase[0], p->phase[1] < 0? '.' : '0' + p->phase[1]);
		}
		fputc('\n', fp);
	}
}
