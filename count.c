#include <assert.h>
#include "hkpriv.h"
#include "kavl.h"
#include "klist.h"

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
	unsigned cl;
	left = 0;
	kavl_insert(nei2, &root, &a[0], 0);
	for (i = 1; i < n_pairs; ++i) {
		for (j = left; j < i; ++j) {
			unsigned cm;
			if (a[i].pos1 - a[j].pos1 < r1) break;
			kavl_erase(nei2, &root, &a[j], &cm);
			a[j].n += count_in_tree2(root, a[j].pos2, r2, &cl);
			assert(cm >= cl + 1);
			a[j].n_corner = cm - 1 - cl;
		}
		left = j;
		assert(i - left == kavl_size(head, root));
		a[i].n = count_in_tree2(root, a[i].pos2, r2, &cl);
		kavl_insert(nei2, &root, &a[i], 0);
	}
	for (j = left; j < n_pairs; ++j) {
		kavl_erase(nei2, &root, &a[j], 0);
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

/******************
 * Cluster by nei *
 ******************/

static void hk_select_by_nei_core(int32_t n_pairs, struct cnt_nei2_aux *a, int radius)
{
	struct cnt_nei2_aux *root = 0;
	int32_t i, j, left, m_del = 0, n_del = 0, *del = 0;
	left = 0;
	kavl_insert(nei2, &root, &a[0], 0);
	for (i = 1; i < n_pairs; ++i) {
		struct cnt_nei2_aux *p, t;
		int to_add = 1;
		for (j = left; j < i; ++j) {
			if (a[i].pos1 - a[j].pos1 < radius) break;
			p = kavl_erase(nei2, &root, &a[j], 0);
			if (p) p->n_corner = 1;
		}
		left = j;
		n_del = 0;
		if (root) {
			kavl_itr_t(nei2) itr;
			t.pos2 = a[i].pos2 > radius? a[i].pos2 - radius : 0;
			t.i = -1;
			kavl_itr_find(nei2, root, &t, &itr);
			while ((p = (struct cnt_nei2_aux*)kavl_at(&itr)) != 0) {
				if (p->pos2 > a[i].pos2 + radius) break;
				if (p->n >= a[i].n) { // there is a better overlapping contact in the tree; then don't add
					to_add = 0;
					break;
				} else {
					if (m_del == n_del)
						EXPAND(del, m_del);
					del[n_del++] = p->i;
				}
				if (!kavl_itr_next(nei2, &itr)) break;
			}
		}
		if (to_add) {
			for (j = 0; j < n_del; ++j)
				kavl_erase(nei2, &root, &a[del[j]], 0);
			kavl_insert(nei2, &root, &a[i], 0);
		}
	}
	while (root) {
		struct cnt_nei2_aux *p;
		p = kavl_erase_first(nei2, &root);
		p->n_corner = 1;
	}
	free(del);
}

int32_t hk_select_by_nei(int32_t n_pairs, struct hk_pair *pairs, int radius, int by_qloop)
{
	int32_t k, st, n_pairs_new = 0;
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
				a[i].n = by_qloop? (int)(p->_.qloop * 100.0 + .499) : p->n_nei;
				a[i].i = i;
			}
			hk_select_by_nei_core(n, a, radius);
			for (i = 0; i < n; ++i)
				if (a[i].n_corner)
					pairs[n_pairs_new++] = pairs[st + a[i].i];
			free(a);
			st = k;
		}
	}
	return n_pairs_new;
}

/***************************
 * Compute expected counts *
 ***************************/

#include "ksort.h"

struct ecnt_aux {
	uint64_t pos;
	int32_t i, n;
};

#define ecnt_key(p) ((p).pos)
KRADIX_SORT_INIT(ecnt, struct ecnt_aux, ecnt_key, 8)

void hk_expected_count(int32_t n_pairs, struct hk_pair *pairs, int r1, int r2)
{
	int32_t i, j, left;
	struct ecnt_aux *a;
	double tmp = 1.0 / n_pairs;
	a = CALLOC(struct ecnt_aux, n_pairs);
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		a[i].i = i;
		a[i].pos = p->chr>>32<<32 | hk_ppos1(p);
	}
	radix_sort_ecnt(a, a + n_pairs);
	for (i = left = 0; i < n_pairs; ++i) {
		struct ecnt_aux *p = &a[i];
		for (j = left; j < i; ++j) {
			if (p->pos - a[j].pos < r1) break;
			a[j].n += i - j;
		}
		left = j;
		p->n = i - left;
	}
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		pairs[a[i].i]._.ecnt.n[0] = a[i].n;
		a[i].i = i;
		a[i].pos = p->chr<<32 | hk_ppos2(p);
		a[i].n = 0;
	}
	radix_sort_ecnt(a, a + n_pairs);
	for (i = left = 0; i < n_pairs; ++i) {
		struct ecnt_aux *p = &a[i];
		for (j = left; j < i; ++j) {
			if (p->pos - a[j].pos < r2) break;
			a[j].n += i - j;
		}
		left = j;
		p->n = i - left;
	}
	for (i = 0; i < n_pairs; ++i)
		pairs[a[i].i]._.ecnt.n[1] = a[i].n;
	free(a);
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		p->_.ecnt.e = tmp * p->_.ecnt.n[0] * p->_.ecnt.n[1];
	}
}
