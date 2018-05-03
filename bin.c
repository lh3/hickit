#include <string.h>
#include <assert.h>
#include "hkpriv.h"
#include "khash.h"
#include "ksort.h"

#define bpair_lt(a, b) ((a).bid[0] < (b).bid[0] || ((a).bid[0] == (b).bid[0] && (a).bid[1] < (b).bid[1]))
KSORT_INIT(bpair, struct hk_bpair, bpair_lt)

void hk_bpair_sort(int32_t n_pairs, struct hk_bpair *pairs)
{
	ks_introsort_bpair(n_pairs, pairs);
}

struct cnt_aux {
	int32_t bid[2];
	float c[4];
};

static inline uint32_t hash_pair(struct cnt_aux x)
{
	return (uint32_t)hash64((uint64_t)x.bid[0] << 32 | x.bid[1]);
}

#define pair_eq(a, b) ((a).bid[0] == (b).bid[0] && (a).bid[1] == (b).bid[1])
KHASH_INIT(bin_cnt, struct cnt_aux, char, 0, hash_pair, pair_eq)

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
}

void hk_bmap_set_offcnt(struct hk_bmap *m)
{
	int32_t i, off = 0;
	assert(m && m->d && m->beads);
	m->offcnt = CALLOC(uint64_t, m->d->n + 1);
	for (i = 1; i <= m->n_beads; ++i) {
		if (i == m->n_beads || m->beads[i].chr != m->beads[i-1].chr) {
			m->offcnt[m->beads[off].chr] = (uint64_t)off << 32 | (i - off);
			off = i;
		}
	}
	m->offcnt[m->d->n] = (uint64_t)off << 32 | 0;
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

struct hk_bmap *hk_bmap_gen(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int size, int min_cnt)
{
	int32_t i, n_del, m_pairs = 0;
	khash_t(bin_cnt) *h;
	struct hk_bmap *m;
	khint_t k;

	m = CALLOC(struct hk_bmap, 1);
	m->ploidy = 1, m->n_full = 0;
	m->d = hk_sd_dup(d, 1, 0);

	hk_bmap_gen_beads_uniform(m, size);
	hk_bmap_set_offcnt(m);

	h = kh_init(bin_cnt);
	for (i = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		struct cnt_aux c;
		int absent;
		khint_t k;
		c.bid[0] = hk_bmap_pos2bid(m, p->chr>>32,      hk_ppos1(p));
		c.bid[1] = hk_bmap_pos2bid(m, (int32_t)p->chr, hk_ppos2(p));
		memcpy(c.c, p->_.p4, sizeof(float) * 4);
		k = kh_put(bin_cnt, h, c, &absent);
		if (!absent) {
			struct cnt_aux *q = &kh_key(h, k);
			q->c[0] += c.c[0], q->c[1] += c.c[1], q->c[2] += c.c[2], q->c[3] += c.c[3];
		}
	}

	for (k = 0, n_del = 0; k < kh_end(h); ++k) {
		struct cnt_aux *q;
		int32_t n;
		if (!kh_exist(h, k)) continue;
		q = &kh_key(h, k);
		n = (int)(q->c[0] + q->c[1] + q->c[2] + q->c[3] + .499f);
		if (n >= min_cnt) {
			int j, max_j = -1;
			float max = -1.0f, scale = 1.0f / n;
			struct hk_bpair *p;
			for (j = 0; j < 4; ++j) {
				q->c[j] *= scale;
				if (max < q->c[j])
					max = q->c[j], max_j = j;
			}
			if (m->n_pairs == m_pairs)
				EXPAND(m->pairs, m_pairs);
			p = &m->pairs[m->n_pairs++];
			p->bid[0] = q->bid[0], p->bid[1] = q->bid[1];
			p->n = n, p->p = max;
			p->phase[0] = max_j>>1&1;
			p->phase[1] = max_j&1;
		}
	}
	kh_destroy(bin_cnt, h);
	hk_bpair_sort(m->n_pairs, m->pairs);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] %d bin pairs\n", __func__, m->n_pairs);
	return m;
}

struct hk_bmap *hk_bmap_dup(const struct hk_bmap *m0, int ploidy, int n_with_homo, int min_cnt, float phase_thres)
{
	struct hk_bmap *m;
	int32_t i, off, chr;
	assert(ploidy >= 1 && ploidy <= 2);
	m = CALLOC(struct hk_bmap, 1);
	m->d = hk_sd_dup(m0->d, ploidy, n_with_homo);
	m->n_beads = (m0->offcnt[n_with_homo]>>32) * ploidy + (m0->n_beads - (m0->offcnt[n_with_homo]>>32));
	m->offcnt = CALLOC(uint64_t, m->d->n + 1);
	m->beads = CALLOC(struct hk_bead, m->n_beads);
	for (i = 0, off = 0, chr = 0; i < m0->d->n; ++i) {
		int32_t off0 = m0->offcnt[i] >> 32;
		int32_t cnt0 = (int32_t)m0->offcnt[i];
		int j, k, pl = i < n_with_homo? ploidy : 1;
		for (j = 0; j < pl; ++j) {
			memcpy(&m->beads[off], &m0->beads[off0], cnt0 * sizeof(struct hk_bead));
			for (k = 0; k < cnt0; ++k) // update chr
				m->beads[off + k].chr = chr;
			m->offcnt[chr++] = (uint64_t)off << 32 | cnt0;
			off += cnt0;
		}
	}
	assert(chr == m->d->n && off == m->n_beads);
	m->offcnt[chr] = (uint64_t)off << 32 | 0;

	m->pairs = CALLOC(struct hk_bpair, m0->n_pairs);
	for (i = 0; i < m0->n_pairs; ++i) {
		struct hk_bpair p = m0->pairs[i];
		int j;
		if (p.n < min_cnt || p.p < phase_thres) continue;
		for (j = 0; j < 2; ++j) {
			int chr_old, chr_new;
			chr_old = m0->beads[p.bid[j]].chr;
			chr_new = chr_old < n_with_homo? ploidy * chr_old + p.phase[j] : n_with_homo * ploidy + (chr_old - n_with_homo);
			p.bid[j] = (m->offcnt[chr_new]>>32) + (p.bid[j] - (m0->offcnt[chr_old]>>32));
		}
		m->pairs[m->n_pairs++] = p;
	}
	hk_bpair_sort(m->n_pairs, m->pairs);
	return m;
}

void hk_bmap_destroy(struct hk_bmap *m)
{
	free(m->x);
	free(m->offcnt);
	free(m->beads);
	free(m->pairs);
	free(m);
}
