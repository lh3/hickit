#include <assert.h>
#include <string.h>
#include "hkpriv.h"
#include "klist.h"
#include "kavl.h"

struct cnt_aux {
	uint64_t end_i;
	uint32_t n_bridges;
	int32_t n_ends; // we can store one "n_bridges+n_ends" instead
	KAVL_HEAD(struct cnt_aux) head;
};

#define cnt_cmp(a, b) (((a)->end_i > (b)->end_i) - ((a)->end_i < (b)->end_i))
KAVL_INIT(cnt, struct cnt_aux, head, cnt_cmp)

#define cnt_free(a)
KMEMPOOL_INIT(cnt, struct cnt_aux, cnt_free)

void hk_pair_count_1chr(int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t i, n, n_ends = 0;
	struct cnt_aux *root = 0;
	kmempool_t(cnt) *mp;
	mp = kmp_init(cnt);
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		int32_t p1 = hk_ppos1(p), p2 = hk_ppos2(p);
		struct cnt_aux *s;
		while (root) {
			struct cnt_aux *t;
			t = kavl_erase_first(cnt, &root);
			if (t->end_i>>32 > p1) { // overlapping _p_
				kavl_insert(cnt, &root, t, 0); // insert back
				break;
			} else {
				struct hk_pair *q = &pairs[(int32_t)t->end_i];
				n = n_ends - t->n_ends - (int32_t)t->n_bridges + 1; // +1 for itself
				assert(n >= 0);
				q->n = n;
				++n_ends;
				kmp_free(cnt, mp, t);
			}
		}
		s = (struct cnt_aux*)kmp_alloc(cnt, mp);
		s->end_i = (uint64_t)p2 << 32 | i;
		s->n_ends = n_ends;
		kavl_insert(cnt, &root, s, &s->n_bridges);
	}
	while (root) {
		struct cnt_aux *t;
		struct hk_pair *q;
		t = kavl_erase_first(cnt, &root);
		q = &pairs[(int32_t)t->end_i];
		n = n_ends - t->n_ends - (int32_t)t->n_bridges + 1;
		assert(n >= 0);
		q->n = n;
		++n_ends;
		kmp_free(cnt, mp, t);
	}
	kmp_destroy(cnt, mp);
}

void hk_pair_count(int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t st, i;
	for (st = 0, i = 1; i <= n_pairs; ++i) {
		if (i == n_pairs || pairs[i].chr != pairs[i-1].chr) {
			if (pairs[st].chr>>32 == (int32_t)pairs[st].chr)
				hk_pair_count_1chr(i - st, &pairs[st]);
			st = i;
		}
	}
} 

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
		f = a[j].mmax_f + ((int32_t)p->n - min_tad_size - area_weight * avg_density * area);
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
		*in_tads_ += pairs[i].n;
		tads[n_tads++] = pairs[i];
	}

	free(a);
	return tads;
}

struct hk_pair *hk_pair2tad(const struct hk_sdict *d, int32_t n_pairs, struct hk_pair *pairs, int min_tad_size, float area_weight, int32_t *n_tads_)
{
	int32_t i, st, *n_tadss, n_tads = 0, tot_in_tads = 0;
	struct hk_pair *tads = 0, **tadss;
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
