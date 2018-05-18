#include <assert.h>
#include "hkpriv.h"
#include "kavl.h"
#include "ksort.h"

struct cnt_nei_aux {
	uint64_t pos1, pos2;
	int32_t i, n;
	KAVL_HEAD(struct cnt_nei_aux) head;
};

#define cnt_nei_key(x) ((x).pos1)
KRADIX_SORT_INIT(nei, struct cnt_nei_aux, cnt_nei_key, 8)

#if 1
#define cnt_nei_cmp(x, y) ((x)->pos2 > (y)->pos2? 1 : (x)->pos2 < (y)->pos2? -1 : (x)->i - (y)->i)
KAVL_INIT(nei, struct cnt_nei_aux, head, cnt_nei_cmp)

static inline int32_t count_in_tree(const struct cnt_nei_aux *root, uint64_t pos, int radius)
{
	struct cnt_nei_aux t;
	unsigned cl, cr;
	if (root == 0) return 0;
	t.i = 0x3fffffff, t.pos2 = pos >> 32 << 32 | ((int32_t)pos > radius? (int32_t)pos - radius : 0);
	kavl_find(nei, root, &t, &cl);
	t.i = -1, t.pos2 = pos + radius;
	kavl_find(nei, root, &t, &cr);
	return cr - cl;
}

static void hk_count_nei_core(int32_t n_pairs, struct cnt_nei_aux *a, int radius)
{
	struct cnt_nei_aux *root = 0;
	int32_t i, j, left;
	left = 0;
	kavl_insert(nei, &root, &a[0], 0);
	for (i = 1; i < n_pairs; ++i) {
		for (j = left; j < i; ++j) {
			if (a[i].pos1 - a[j].pos1 < radius) break;
			kavl_erase(nei, &root, &a[j]);
			a[j].n += count_in_tree(root, a[j].pos2, radius);
		}
		left = j;
		assert(i - left == kavl_size(head, root));
		a[i].n = count_in_tree(root, a[i].pos2, radius);
		kavl_insert(nei, &root, &a[i], 0);
	}
	for (j = left; j < n_pairs; ++j) {
		kavl_erase(nei, &root, &a[j]);
		a[j].n += count_in_tree(root, a[j].pos2, radius);
	}
}
#else
static void hk_count_nei_core(int32_t n_pairs, struct cnt_nei_aux *a, int radius) // simpler but much slower algorithm
{
	int32_t i, j, left;
	left = 0;
	for (i = 1; i < n_pairs; ++i) {
		uint64_t lo, hi;
		for (j = left; j < i; ++j)
			if (a[i].pos1 - a[j].pos1 < radius)
				break;
		left = j;
		lo = a[i].pos2 >> 32 << 32 | ((int32_t)a[i].pos2 > radius? a[i].pos2 - radius : 0);
		hi = a[i].pos2 + radius;
		for (j = left; j < i; ++j)
			if (a[j].pos2 > lo && a[j].pos2 < hi)
				++a[i].n, ++a[j].n;
	}
}
#endif

void hk_pair_count_nei(int32_t n_pairs, struct hk_pair *pairs, int radius)
{
	struct cnt_nei_aux *a;
	int32_t i;
	a = CALLOC(struct cnt_nei_aux, n_pairs);
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		a[i].pos1 = p->chr >> 32 << 32 | p->pos >> 32;
		a[i].pos2 = p->chr << 32 | p->pos << 32 >> 32;
		a[i].i = i;
	}
	radix_sort_nei(a, a + n_pairs);
	hk_count_nei_core(n_pairs, a, radius);
	for (i = 0; i < n_pairs; ++i)
		pairs[a[i].i].n = a[i].n;
	free(a);
}

void hk_bmap_count_nei(struct hk_bmap *m, int radius)
{
	struct cnt_nei_aux *a;
	int32_t i;
	a = CALLOC(struct cnt_nei_aux, m->n_pairs);
	for (i = 0; i < m->n_pairs; ++i) {
		struct hk_bpair *p = &m->pairs[i];
		a[i].pos1 = (uint64_t)m->beads[p->bid[0]].chr << 32 | m->beads[p->bid[0]].st;
		a[i].pos2 = (uint64_t)m->beads[p->bid[1]].chr << 32 | m->beads[p->bid[1]].st;
	}
	hk_count_nei_core(m->n_pairs, a, radius);
	for (i = 0; i < m->n_pairs; ++i)
		m->pairs[a[i].i].n_nei = a[i].n;
	free(a);
}
