#include <string.h>
#include <assert.h>
#include "hkpriv.h"
#include "khash.h"
#include "ksort.h"

#define bpair_lt(a, b) ((a).chr < (b).chr || ((a).chr == (b).chr && (a).pos < (b).pos))
KSORT_INIT(bpair, struct hk_bpair, bpair_lt)

void hk_bpair_sort(int32_t n_pairs, struct hk_bpair *pairs)
{
	ks_introsort_bpair(n_pairs, pairs);
}

static inline uint64_t hash64(uint64_t key)
{
	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return key;
}

struct cnt_aux {
	uint64_t pos[2];
	float c[4];
};

static inline uint32_t hash_pair(struct cnt_aux x)
{
	return (uint32_t)(hash64(x.pos[0]) + hash64(x.pos[1]));
}

#define pair_eq(a, b) ((a).pos[0] == (b).pos[0] && (a).pos[1] == (b).pos[1])
KHASH_INIT(bin_cnt, struct cnt_aux, char, 0, hash_pair, pair_eq)

struct hk_bmap *hk_bmap_gen2(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int size, int n_full, int min_cnt, float phase_thres)
{
	int32_t i, n_del, m_pairs = 0;
	khash_t(bin_cnt) *h;
	struct hk_bmap *m;
	khint_t k;

	m = CALLOC(struct hk_bmap, 1);
	m->d = hk_sd_dup(d, 2, n_full);
	m->ploidy = 2, m->n_full = n_full;

	h = kh_init(bin_cnt);
	for (i = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		struct cnt_aux c;
		int absent;
		khint_t k;
		c.pos[0] = p->chr>>32<<32 | hk_ppos1(p) / size * size;
		c.pos[1] = p->chr<<32     | hk_ppos2(p) / size * size;
		memcpy(c.c, p->_.p4, sizeof(float) * 4);
		k = kh_put(bin_cnt, h, c, &absent);
		if (!absent) {
			struct cnt_aux *q = &kh_key(h, k);
			q->c[0] += c.c[0], q->c[1] += c.c[1], q->c[2] += c.c[2], q->c[3] += c.c[3];
		}
	}

	for (k = 0, n_del = 0; k < kh_end(h); ++k) {
		struct cnt_aux *q;
		int n;
		if (!kh_exist(h, k)) continue;
		q = &kh_key(h, k);
		n = (int)(q->c[0] + q->c[1] + q->c[2] + q->c[3] + .499f);
		if (n >= min_cnt) {
			int32_t j, max_j = -1;
			float max = -1.0f, scale = 1.0f / n;
			for (j = 0; j < 4; ++j) {
				q->c[j] *= scale;
				if (max < q->c[j])
					max = q->c[j], max_j = j;
			}
			if (max >= phase_thres) {
				struct hk_bpair *p;
				int32_t chr[2], pos[2], end[2];
				if (m->n_pairs == m_pairs)
					EXPAND(m->pairs, m_pairs);
				chr[0] = q->pos[0]>>32, chr[1] = q->pos[1]>>32;
				pos[0] = (int32_t)q->pos[0], pos[1] = (int32_t)q->pos[1];
				end[0] = pos[0] + size < d->len[chr[0]]? pos[0] + size : d->len[chr[0]];
				end[1] = pos[1] + size < d->len[chr[1]]? pos[1] + size : d->len[chr[1]];
				chr[0] = chr[0] < n_full? chr[0] * 2 + (max_j>>1&1) : n_full * 2 + (chr[0] - n_full);
				chr[1] = chr[1] < n_full? chr[1] * 2 + (max_j&1) : n_full * 2 + (chr[1] - n_full);
				p = &m->pairs[m->n_pairs++];
				p->chr = (uint64_t)chr[0]<<32 | chr[1];
				p->pos = (uint64_t)pos[0]<<32 | pos[1];
				p->end = (uint64_t)end[0]<<32 | end[1];
				p->n = n, p->p = max;
			}
		}
	}
	kh_destroy(bin_cnt, h);
	hk_bpair_sort(m->n_pairs, m->pairs);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] %d bin pairs\n", __func__, m->n_pairs);
	return m;
}

void hk_bmap_destroy(struct hk_bmap *m)
{
	free(m->pairs);
	free(m);
}
