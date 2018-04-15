#include <assert.h>
#include "hkpriv.h"
#include "kavl.h"

struct cnt_aux {
	uint64_t end_i;
	uint32_t n_bridges;
	int32_t n_ends; // we can store one "n_bridges+n_ends" instead
	KAVL_HEAD(struct cnt_aux) head;
};

#define cnt_cmp(a, b) (((a)->end_i > (b)->end_i) - ((a)->end_i < (b)->end_i))
KAVL_INIT(cnt, struct cnt_aux, head, cnt_cmp)

void hk_pair_count_1chr(int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t i, n_ends = 0;
	struct cnt_aux *root = 0;
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
				q->n = n_ends - t->n_ends - (int32_t)t->n_bridges;
				assert(q->n >= 0);
				++n_ends;
				free(t);
			}
		}
		s = CALLOC(struct cnt_aux, 1);
		s->end_i = (uint64_t)p2 << 32 | i;
		s->n_ends = n_ends;
		kavl_insert(cnt, &root, s, &s->n_bridges);
	}
	while (root) {
		struct cnt_aux *t;
		struct hk_pair *q;
		t = kavl_erase_first(cnt, &root);
		q = &pairs[(int32_t)t->end_i];
		q->n = n_ends - t->n_ends - (int32_t)t->n_bridges;
		++n_ends;
		free(t);
	}
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
