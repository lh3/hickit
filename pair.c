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
					p->n = 0, p->offset = -1;
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

int32_t hk_pair_filter(int n_pairs, struct hk_pair *pairs, int min_dist)
{
	int32_t i, n;
	for (i = n = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (hk_intra(p) && hk_ppos2(p) - hk_ppos1(p) < min_dist)
			continue;
		pairs[n++] = pairs[i];
	}
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] filtered %d out of %d pairs\n", __func__, n_pairs - n, n_pairs);
	return n;
}

int32_t hk_mask_by_tad(int32_t n_tads, const struct hk_pair *tads, int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t i, j, k;
	for (i = j = k = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		int kept = 1;
		if (p->chr>>32 == (int32_t)p->chr) { // intra-chromosomal pairs
			const struct hk_pair *t;
			int32_t p1 = hk_ppos1(p);
			while (j < n_tads && (tads[j].chr < p->chr || (tads[j].chr == p->chr && hk_ppos2(&tads[j]) <= p1)))
				++j;
			if (j == n_tads) break;
			t = &tads[j];
			if (p->chr == t->chr && hk_ppos1(t) <= p1 && hk_ppos2(p) <= hk_ppos2(t)) // contained in TAD
				kept = 0;
		}
		if (kept) pairs[k++] = pairs[i];
	}
	for (; i < n_pairs; ++i) // copy the rest over
		pairs[k++] = pairs[i];
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

void hk_print_pair(FILE *fp, const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs)
{
	int32_t i;
	fprintf(fp, "## pairs format v1.0\n");
	fprintf(fp, "#sorted: chr1-chr2-pos1-pos2\n");
	fprintf(fp, "#shape: upper triangle\n");
	hk_print_chr(fp, d);
	fprintf(fp, "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n");
	for (i = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		fprintf(fp, ".\t%s\t%d\t%s\t%d\t%c\t%c\t%d\n", d->name[p->chr>>32], (int32_t)(p->pos>>32),
				d->name[(int32_t)p->chr], (int32_t)p->pos,
				p->strand[0] > 0? '+' : p->strand[0] < 0? '-' : '.',
				p->strand[1] > 0? '+' : p->strand[1] < 0? '-' : '.', p->n);
	}
}

/***************
 * TAD calling *
 ***************/

struct dp_aux {
	float f;
	int32_t j;
};

static struct hk_pair *hk_tad_call1(int32_t n_pairs, struct hk_pair *pairs, int max_radius, float area_weight, int min_back, int32_t *n_tads_, int *in_band, int *in_tads)
{
	int32_t i, n_a, max_max_i, n_tads, min = INT32_MAX, max = INT32_MIN;
	struct hk_pair *a, *tads;
	float area_tot, aa, max_max_f;
	struct dp_aux *u;

	*in_band = *in_tads = 0;

	// generate a[]
	for (i = n_a = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (hk_ppos2(p) - hk_ppos1(p) <= max_radius)
			++n_a;
	}
	a = MALLOC(struct hk_pair, n_a);
	for (i = n_a = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		int32_t p1 = hk_ppos1(p), p2 = hk_ppos2(p);
		if (p2 - p1 <= max_radius) {
			a[n_a] = *p;
			a[n_a++].n = 0;
			min = min < p1? min : p1;
			max = max > p2? max : p2;
		}
	}
	area_tot = 1e-12 * (max - min + 0.5 * max_radius) * max_radius; // #pairs per square-million
	aa = area_weight * n_a / area_tot;

	// count
	for (i = 1; i < n_a; ++i) {
		struct hk_pair *q = &a[i];
		int32_t j, q1 = hk_ppos1(q), q2 = hk_ppos2(q);
		for (j = i - 1; j >= 0; --j) {
			struct hk_pair *p = &a[j];
			int32_t p1 = hk_ppos1(p), p2 = hk_ppos2(p);
			if (q2 - p1 > max_radius) break;
			if (q1 != p1 && q2 < p2) ++p->n;
		}
	}

	// DP
	u = CALLOC(struct dp_aux, n_a);
	max_max_f = -2e30f, max_max_i = -i;
	for (i = n_a - 1; i >= 0; --i) {
		struct hk_pair *p = &a[i];
		int32_t j, max_j = n_a, p1 = hk_ppos1(p), p2 = hk_ppos2(p);
		float area = 0.5e-12f * (p2 - p1) * (p2 - p1);
		//float max_f = p->n - aa * (area + 0.5e-12f * (max - p2) * (max - p2));
		float max_f = p->n - aa * area;
		for (j = i + 1; j < n_a; ++j) {
			struct hk_pair *q = &a[j];
			int32_t q1 = hk_ppos1(q), q2;
			float f;
			if (q1 - p1 > max_radius && j - i >= min_back) break;
			if (p2 >= q1) continue;
			q2 = hk_ppos2(q);
			//f = u[j].f + p->n - aa * (area + 0.5e-12f * (q1 - p2) * (q1 - p2));
			f = u[j].f + p->n - aa * area;
			if (f > max_f)
				max_f = f, max_j = j;
		}
		u[i].f = max_f, u[i].j = max_j;
		if (max_f > max_max_f)
			max_max_f = max_f, max_max_i = i;
	}

	for (i = max_max_i, n_tads = 0; i < n_a; i = u[i].j)
		++n_tads;
	*n_tads_ = n_tads;
	tads = CALLOC(struct hk_pair, n_tads);
	for (i = max_max_i, n_tads = 0; i < n_a; i = u[i].j) {
		*in_tads += a[i].n;
		tads[n_tads++] = a[i];
	}
	*in_band = n_a;

	free(u);
	free(a);
	return tads;
}

struct hk_pair *hk_pair2tad_slow(const struct hk_sdict *d, int32_t n_pairs, struct hk_pair *pairs, int max_radius, float area_weight, int32_t *n_tads_)
{
	int32_t i, st, *n_tadss, n_tads = 0, tot_in_band = 0, tot_in_tads = 0;
	struct hk_pair *tads = 0, **tadss;
	tadss = CALLOC(struct hk_pair*, d->n);
	n_tadss = CALLOC(int32_t, d->n);
	for (st = 0, i = 1; i <= n_pairs; ++i) {
		if (i == n_pairs || pairs[i].chr != pairs[i-1].chr) {
			if (pairs[st].chr>>32 == (int32_t)pairs[st].chr) {
				int32_t chr = (int32_t)pairs[st].chr, in_band, in_tads;
				tadss[chr] = hk_tad_call1(i - st, &pairs[st], max_radius, area_weight, 64, &n_tadss[chr], &in_band, &in_tads);
				n_tads += n_tadss[chr];
				tot_in_band += in_band;
				tot_in_tads += in_tads;
			}
			st = i;
		}
	}
	tads = MALLOC(struct hk_pair, n_tads);
	for (i = 0, n_tads = 0; i < d->n; ++i) {
		if (n_tadss[i] > 0) {
			memcpy(&tads[n_tads], tadss[i], n_tadss[i] * sizeof(struct hk_pair));
			n_tads += n_tadss[i];
			free(tadss[i]);
		}
	}
	free(tadss);
	free(n_tadss);
	*n_tads_ = n_tads;
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] %d pairs in band, of which %d in %d TADs (%.2f%%)\n", __func__,
				tot_in_band, tot_in_tads, n_tads, 100.0 * tot_in_tads / tot_in_band);
	return tads;
}
