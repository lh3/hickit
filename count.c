#include <assert.h>
#include "hkpriv.h"
#include "kavl.h"
#include "klist.h"
#include "ksort.h"

/*************************
 * Count contained pairs *
 *************************/

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
				q->n_ctn = n;
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
		q->n_ctn = n;
		++n_ends;
		kmp_free(cnt, mp, t);
	}
	kmp_destroy(cnt, mp);
}

void hk_pair_count_contained(int32_t n_pairs, struct hk_pair *pairs)
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

/***************************
 * Count neighboring pairs *
 ***************************/

struct cnt_nei2_aux {
	int32_t pos1, pos2; // the algorithm can be modified to avoid pos1, which will save 8 bytes
	int32_t i, n, n_corner;
	KAVL_HEAD(struct cnt_nei2_aux) head;
};

#define cnt_nei2_cmp(x, y) ((x)->pos2 > (y)->pos2? 1 : (x)->pos2 < (y)->pos2? -1 : (x)->i - (y)->i)
KAVL_INIT(nei2, struct cnt_nei2_aux, head, cnt_nei2_cmp)

static inline int32_t count_in_tree2(const struct cnt_nei2_aux *root, int32_t pos, int radius, unsigned *cl)
{
	struct cnt_nei2_aux t;
	unsigned cr;
	*cl = 0;
	if (root == 0) return 0;
	t.i = 0x3fffffff, t.pos2 = pos > radius? pos - radius : 0;
	kavl_find(nei2, root, &t, cl);
	t.i = -1, t.pos2 = pos + radius;
	kavl_find(nei2, root, &t, &cr);
	return cr - *cl;
}

static void hk_count_nei2_core(int32_t n_pairs, struct cnt_nei2_aux *a, int r1, int r2)
{
	struct cnt_nei2_aux *root = 0;
	int32_t i, j, left;
	unsigned cl, cm;
	left = 0;
	kavl_insert(nei2, &root, &a[0], 0);
	for (i = 1; i < n_pairs; ++i) {
		for (j = left; j < i; ++j) {
			if (a[i].pos1 - a[j].pos1 < r1) break;
			kavl_erase(nei2, &root, &a[j]);
			a[j].n += count_in_tree2(root, a[j].pos2, r2, &cl);
		}
		left = j;
		assert(i - left == kavl_size(head, root));
		a[i].n = count_in_tree2(root, a[i].pos2, r2, &cl);
		kavl_insert(nei2, &root, &a[i], &cm);
		assert(cm >= cl);
		a[i].n_corner = cm - cl;
	}
	for (j = left; j < n_pairs; ++j) {
		kavl_erase(nei2, &root, &a[j]);
		a[j].n += count_in_tree2(root, a[j].pos2, r2, &cl);
	}
}

void hk_pair_count_nei(int32_t n_pairs, struct hk_pair *pairs, int r1, int r2)
{
	int32_t k, st;
	for (k = 1; k < n_pairs; ++k)
		if (pairs[k-1].chr > pairs[k].chr || (pairs[k-1].chr == pairs[k].chr && pairs[k-1].pos > pairs[k].pos))
			break;
	assert(k == n_pairs); // otherwise, pairs[] is not sorted
	for (k = 1, st = 0; k <= n_pairs; ++k) {
		if (k == n_pairs || pairs[k-1].chr != pairs[k].chr) {
			struct cnt_nei2_aux *a;
			int i, n = k - st;
			a = CALLOC(struct cnt_nei2_aux, n);
			for (i = 0; i < n; ++i) {
				struct hk_pair *p = &pairs[st + i];
				a[i].pos1 = hk_ppos1(p);
				a[i].pos2 = hk_ppos2(p);
				a[i].i = i;
			}
			hk_count_nei2_core(n, a, r1, r2);
			for (i = 0; i < n; ++i) {
				pairs[st + a[i].i].n_nei = a[i].n;
				pairs[st + a[i].i].n_nei_corner = a[i].n_corner;
			}
			free(a);
			st = k;
		}
	}
}
