#include <string.h>
#include <assert.h>
#include "hkpriv.h"
#include "khash.h"
#include "ksort.h"

#define bpair_lt(a, b) ((a).bid[0] < (b).bid[0] || ((a).bid[0] == (b).bid[0] && (a).bid[1] < (b).bid[1]))
KSORT_INIT(bpair, struct hk_bpair, bpair_lt)

KSORT_INIT_GENERIC(int32_t)

void hk_bpair_sort(int32_t n_pairs, struct hk_bpair *pairs)
{
	ks_introsort_bpair(n_pairs, pairs);
}

struct cnt_aux {
	int32_t bid[2];
	int32_t n;
	float p_sum;
};

static inline uint32_t hash_pair(struct cnt_aux x)
{
	return (uint32_t)hash64((uint64_t)x.bid[0] << 32 | x.bid[1]);
}

#define pair_eq(a, b) ((a).bid[0] == (b).bid[0] && (a).bid[1] == (b).bid[1])
KHASH_INIT(bin_cnt, struct cnt_aux, char, 0, hash_pair, pair_eq)

void hk_bmap_set_offcnt(struct hk_bmap *m)
{
	int32_t i, off = 0;
	assert(m && m->d && m->beads);
	if (m->offcnt) free(m->offcnt);
	m->offcnt = CALLOC(uint64_t, m->d->n + 1);
	for (i = 1; i <= m->n_beads; ++i) {
		if (i == m->n_beads || m->beads[i].chr != m->beads[i-1].chr) {
			m->offcnt[m->beads[off].chr] = (uint64_t)off << 32 | (i - off);
			off = i;
		}
	}
	m->offcnt[m->d->n] = (uint64_t)off << 32 | 0;
}

void hk_bmap_gen_beads_uniform(struct hk_bmap *m, int size)
{
	int32_t i, m_beads = 0;
	assert(m && m->d);
	for (i = 0; i < m->d->n; ++i) {
		int32_t st = 0, len = m->d->len[i];
		while (st < len) {
			int32_t l = size * 1.5 < len - st? size : len - st;
			struct hk_bead *p;
			if (m->n_beads == m_beads)
				EXPAND(m->beads, m_beads);
			p = &m->beads[m->n_beads++];
			p->chr = i, p->st = st, p->en = st + l;
			st += l;
		}
	}
	hk_bmap_set_offcnt(m);
}

int hk_bmap_pos2bid(const struct hk_bmap *m, int32_t chr, int32_t pos)
{
	int32_t lo, hi;
	lo = m->offcnt[chr] >> 32;
	hi = lo + (int32_t)m->offcnt[chr] - 1;
	assert(pos < m->beads[hi].en);
	while (lo <= hi) {
		int32_t mid = (lo + hi) / 2;
		const struct hk_bead *p = &m->beads[mid];
		if (p->st <= pos && pos < p->en) return mid;
		else if (p->st < pos) lo = mid + 1;
		else hi = mid - 1;
	}
	abort(); // if here, it is a bug
	return -1;
}

static khash_t(bin_cnt) *hash_count(struct hk_bmap *m, int32_t n_pairs, const struct hk_pair *pairs)
{
	int32_t i;
	khash_t(bin_cnt) *h;
	h = kh_init(bin_cnt);
	for (i = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		struct cnt_aux c;
		int absent;
		khint_t k;
		c.bid[0] = hk_bmap_pos2bid(m, p->chr>>32,      hk_ppos1(p));
		c.bid[1] = hk_bmap_pos2bid(m, (int32_t)p->chr, hk_ppos2(p));
		c.n = 1, c.p_sum = p->_.phased_prob;
		k = kh_put(bin_cnt, h, c, &absent);
		if (!absent) {
			struct cnt_aux *q = &kh_key(h, k);
			++q->n, q->p_sum += c.p_sum;
		}
	}
	return h;
}

void hk_bmap_merge_beads(struct hk_bmap *m, int32_t n_pairs, const struct hk_pair *pairs)
{
	int32_t i, n, *cnt;
	khash_t(bin_cnt) *h;
	khint_t k;

	h = hash_count(m, n_pairs, pairs);
	cnt = CALLOC(int32_t, m->n_beads);
	for (k = 0; k < kh_end(h); ++k) {
		struct cnt_aux *q;
		if (!kh_exist(h, k)) continue;
		q = &kh_key(h, k);
		if (q->bid[0] == q->bid[1]) continue;
		if (m->beads[q->bid[0]].chr == m->beads[q->bid[1]].chr && q->bid[1] == q->bid[0] + 1) continue;
		++cnt[q->bid[0]];
		++cnt[q->bid[1]];
	}

	for (i = n = 0; i < m->n_beads; ++i) {
		if (cnt[i] < 1) {
			if (i == 0 || m->beads[i-1].chr != m->beads[i].chr)
				m->beads[n++] = m->beads[i];
			else m->beads[n-1].en = m->beads[i].en;
		} else m->beads[n++] = m->beads[i];
	}
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] %d => %d\n", __func__, m->n_beads, n);
	m->n_beads = n;

	free(cnt);
	kh_destroy(bin_cnt, h);
	hk_bmap_set_offcnt(m);
}

struct hk_bmap *hk_bmap_gen(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int size)
{
	int32_t n_del, m_pairs = 0;
	khash_t(bin_cnt) *h;
	struct hk_bmap *m;
	khint_t k;

	m = CALLOC(struct hk_bmap, 1);
	m->d = hk_sd_dup(d, 1, 0);

	hk_bmap_gen_beads_uniform(m, size);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] generated %d beads\n", __func__, m->n_beads);
	hk_bmap_merge_beads(m, n_pairs, pairs);
	hk_bmap_merge_beads(m, n_pairs, pairs);

	h = hash_count(m, n_pairs, pairs);
	for (k = 0, n_del = 0; k < kh_end(h); ++k) {
		struct hk_bpair *p;
		struct cnt_aux *q;
		if (!kh_exist(h, k)) continue;
		q = &kh_key(h, k);
		if (m->n_pairs == m_pairs)
			EXPAND(m->pairs, m_pairs);
		p = &m->pairs[m->n_pairs++];
		p->bid[0] = q->bid[0], p->bid[1] = q->bid[1];
		p->n = q->n, p->p = q->p_sum / q->n;
	}
	kh_destroy(bin_cnt, h);
	hk_bpair_sort(m->n_pairs, m->pairs);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] %d bin pairs\n", __func__, m->n_pairs);
	return m;
}

void hk_bmap_destroy(struct hk_bmap *m)
{
	if (m->d) hk_sd_destroy(m->d);
	free(m->x);
	free(m->offcnt);
	free(m->beads);
	free(m->pairs);
	free(m);
}
