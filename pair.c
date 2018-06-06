#include <assert.h>
#include <string.h>
#include "hkpriv.h"
#include "ksort.h"
#define pair_lt(a, b) ((a).chr < (b).chr || ((a).chr == (b).chr && (a).pos < (b).pos))
KSORT_INIT(pair, struct hk_pair, pair_lt)

int hk_verbose = 3;

void hk_popt_init(struct hk_popt *c)
{
	memset(c, 0, sizeof(struct hk_popt));
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

/***************************************
 * Assorted utilities to process pairs *
 ***************************************/

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
					memset(p, 0, sizeof(struct hk_pair));
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
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] generated %d pairs\n", __func__, n_pairs);
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

void hk_mark_by_tad(int32_t n_tads, const struct hk_pair *tads, int32_t n_pairs, struct hk_pair *pairs)
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
		pairs[i].tad_marked = !kept;
		if (!kept) ++n_masked;
	}
	free(a);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] masked %d out of %d pairs\n", __func__, n_masked, n_pairs);
}

struct hk_map *hk_pair_split_phase(const struct hk_map *m, float phase_thres)
{
	int32_t i, m_ppairs = 0, *ploidy_XY, *old2new;
	struct hk_map *p;

	ploidy_XY = hk_sd_ploidy_XY(m->d, 0);
	old2new = CALLOC(int32_t, m->d->n);
	for (i = 1; i < m->d->n; ++i)
		old2new[i] = old2new[i-1] + (ploidy_XY[i-1]>>8);

	p = CALLOC(struct hk_map, 1);
	p->d = hk_sd_split_phase(m->d, ploidy_XY);
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

int32_t hk_pair_filter_isolated(int32_t n_pairs, struct hk_pair *pairs, int32_t max_radius, int32_t min_cnt, float drop_frac)
{
	int32_t i, k, *cnt, min;
	hk_pair_count_nei(n_pairs, pairs, max_radius);
	cnt = CALLOC(int32_t, n_pairs);
	for (i = 0; i < n_pairs; ++i)
		cnt[i] = pairs[i].n_nei;
	min = ks_ksmall_int32_t(n_pairs, cnt, n_pairs * drop_frac);
	min = min > min_cnt? min : min_cnt;
	free(cnt);
	for (i = k = 0; i < n_pairs; ++i)
		if (pairs[i].n_nei >= min)
			pairs[k++] = pairs[i];
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] threshold: %d; filtered out %d/%d isolated pairs\n", __func__, min, n_pairs - k, n_pairs);
	return k;
}

int32_t hk_pair_filter_close_legs(int32_t n_pairs, struct hk_pair *pairs, int min_dist)
{
	int32_t i, k;
	for (i = k = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (!hk_intra(p) || hk_ppos2(p) - hk_ppos1(p) >= min_dist || p->strand[0] * p->strand[1] < 0)
			pairs[k++] = *p;
	}
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] filtered out %d/%d pairs with close legs\n", __func__, n_pairs - k, n_pairs);
	return k;
}

/***************
 * TAD calling *
 ***************/

struct tad_aux {
	float mmax_f;
	int32_t i, mmax_i;
};

static struct hk_pair *hk_tad_call1(int32_t n_pairs, struct hk_pair *pairs, int min_tad_size, float area_weight, int32_t *n_tads_, int32_t *in_tads_)
{
	int32_t i, min_pos, max_pos, mmax_i, n_tads;
	float mmax_f, avg_density;
	struct tad_aux *a;
	struct hk_pair *tads;

	min_pos = INT32_MAX, max_pos = INT32_MIN;
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		int32_t p1 = hk_ppos1(p), p2 = hk_ppos2(p);
		min_pos = min_pos < p1? min_pos : p1;
		max_pos = max_pos > p2? max_pos : p2;
	}
	avg_density = n_pairs / (0.5e-12f * (max_pos - min_pos) * (max_pos - min_pos));

	a = CALLOC(struct tad_aux, n_pairs + 1);
	a[n_pairs].mmax_f = mmax_f = 0.0f;
	a[n_pairs].mmax_i = mmax_i = n_pairs;
	for (i = n_pairs - 1; i >= 0; --i) {
		struct hk_pair *p = &pairs[i];
		int32_t p1 = hk_ppos1(p), p2 = hk_ppos2(p);
		float f, area = 0.5e-12f * (p2 - p1) * (p2 - p1);
		int32_t j = n_pairs, lo = i + 1, hi = n_pairs - 1;
		while (lo <= hi) { // binary search for the nearest pair that starts at or after _p2_
			int32_t mid = (lo + hi) / 2;
			if (hk_ppos1(&pairs[mid]) < p2) {
				lo = mid + 1;
				if (lo >= n_pairs || hk_ppos1(&pairs[lo]) >= p2) {
					j = mid;
					break;
				}
			} else {
				hi = mid - 1;
				if (hi < i + 1 || hk_ppos1(&pairs[hi]) < p2) {
					j = mid;
					break;
				}
			}
		}
		a[i].i = a[j].mmax_i;
		f = a[j].mmax_f + ((int32_t)p->n_ctn - min_tad_size - area_weight * avg_density * area);
		if (f >= mmax_f)
			mmax_f = f, mmax_i = i;
		a[i].mmax_f = mmax_f;
		a[i].mmax_i = mmax_i;
	}

	for (i = mmax_i, n_tads = 0; i < n_pairs; i = a[i].i)
		++n_tads;
	*n_tads_ = n_tads;
	tads = CALLOC(struct hk_pair, n_tads);
	for (i = mmax_i, n_tads = 0i, *in_tads_ = 0; i < n_pairs; i = a[i].i) {
		*in_tads_ += pairs[i].n_ctn;
		tads[n_tads++] = pairs[i];
	}

	free(a);
	return tads;
}

struct hk_pair *hk_pair2tad(const struct hk_sdict *d, int32_t n_pairs, struct hk_pair *pairs, int min_tad_size, float area_weight, int32_t *n_tads_)
{
	int32_t i, st, *n_tadss, n_tads = 0, tot_in_tads = 0;
	struct hk_pair *tads = 0, **tadss;
	for (i = 0; i < n_pairs; ++i)
		if (pairs[i].n_ctn > 0) break;
	if (i == n_pairs) hk_pair_count_contained(n_pairs, pairs);
	tadss = CALLOC(struct hk_pair*, d->n);
	n_tadss = CALLOC(int32_t, d->n);
	for (st = 0, i = 1; i <= n_pairs; ++i) {
		if (i == n_pairs || pairs[i].chr != pairs[i-1].chr) {
			if (pairs[st].chr>>32 == (int32_t)pairs[st].chr) {
				int32_t chr = (int32_t)pairs[st].chr, in_tads;
				tadss[chr] = hk_tad_call1(i - st, &pairs[st], min_tad_size, area_weight, &n_tadss[chr], &in_tads);
				n_tads += n_tadss[chr];
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
		fprintf(stderr, "[M::%s] %d pairs (%.2f%%) in %d TADs\n", __func__,
				tot_in_tads, 100.0 * tot_in_tads / n_pairs, n_tads);
	return tads;
}
