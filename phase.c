#include <math.h>
#include "hkpriv.h"

/*
float hk_pair_weight(const struct hk_pair *a, const struct hk_pair *b, int32_t max, float alpha, float beta)
{
	int32_t y, z;
	float yf, zf, d;
	if (a->chr != b->chr) return 0.0f;
	y = hk_ppos1(b) > hk_ppos1(a)? hk_ppos1(b) - hk_ppos1(a) : hk_ppos1(a) - hk_ppos1(b);
	z = hk_ppos2(b) > hk_ppos2(a)? hk_ppos2(b) - hk_ppos2(a) : hk_ppos2(a) - hk_ppos2(b);
	if (y > max || z > max)
		return 0.0f;
	yf = y / (float)max;
	zf = z / (float)max;
	d = 0.5f * (yf * yf + zf * zf + (alpha + alpha - 1.0f) * (yf - zf) * (yf - zf));
	return d >= 1.0f? 0.0f : expf(-beta * d);
}
*/

static float dist2weight(int32_t d, int32_t max, float beta)
{
	float x = (float)d / max;
	return expf(-beta * x * x);
}

void hk_nei_weight(struct hk_nei *n, int32_t max_radius, float beta)
{
	int32_t i;
	for (i = 0; i < n->n_pairs; ++i) {
		int64_t off = n->offcnt[i] >> 16;
		int32_t cnt = n->offcnt[i] & 0xffff, j;
		for (j = 0; j < cnt; ++j) {
			struct hk_nei1 *n1 = &n->nei[off + j];
			n1->_.w = dist2weight(n1->_.d, max_radius, beta);
		}
	}
}

struct phase_aux {
	float phase[2];
};

void hk_nei_phase(struct hk_nei *n, struct hk_pair *pairs, int n_iter, float pseudo_cnt)
{
	struct phase_aux *a[2], *cur, *pre, *tmp;
	int32_t i, iter;

	a[0] = CALLOC(struct phase_aux, n->n_pairs * 2);
	a[1] = a[0] + n->n_pairs;
	cur = a[0], pre = a[1];

	// initialize
	for (i = 0; i < n->n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		cur[i].phase[0] = p->phase[0] < 0? 0.5f : p->phase[0];
		cur[i].phase[1] = p->phase[1] < 0? 0.5f : p->phase[1];
	}
	tmp = cur, cur = pre, pre = tmp;

	// EM, I think
	for (iter = 0; iter < n_iter; ++iter) {
		for (i = 0; i < n->n_pairs; ++i) {
			struct hk_pair *p = &pairs[i];
			int64_t off = n->offcnt[i] >> 16;
			int32_t cnt = n->offcnt[i] & 0xffff, j;
			double c[2], s;
			c[0] = c[1] = 0.5 * pseudo_cnt, s = pseudo_cnt;
			for (j = 0; j < cnt; ++j) {
				struct hk_nei1 *n1 = &n->nei[off + j];
				c[0] += n1->_.w * pre[n1->i].phase[0];
				c[1] += n1->_.w * pre[n1->i].phase[1];
				s += n1->_.w;
			}
			cur[i].phase[0] = p->phase[0] >= 0? p->phase[0] : c[0] / s;
			cur[i].phase[1] = p->phase[1] >= 0? p->phase[1] : c[1] / s;
		}
		tmp = cur, cur = pre, pre = tmp;
		fprintf(stderr, "%d\n", iter);
	}

	// write back
	for (i = 0; i < n->n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		p->phase[0] = (int)(pre[i].phase[0] * 100.0f + .499f);
		p->phase[1] = (int)(pre[i].phase[1] * 100.0f + .499f);
	}
	free(a[0]);
}

/***************************
 * Random number generator *
 ***************************/

typedef struct {
	uint64_t s[2];
	double n_gset;
	int n_iset;
	volatile int lock;
} kad_rng_t;

static kad_rng_t kad_rng_dat = { {0x50f5647d2380309dULL, 0x91ffa96fc4c62cceULL}, 0.0, 0, 0 };

static inline uint64_t kad_splitmix64(uint64_t x)
{
	uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static inline uint64_t kad_xoroshiro128plus_next(kad_rng_t *r)
{
	const uint64_t s0 = r->s[0];
	uint64_t s1 = r->s[1];
	const uint64_t result = s0 + s1;
	s1 ^= s0;
	r->s[0] = (s0 << 55 | s0 >> 9) ^ s1 ^ (s1 << 14);
	r->s[1] = s0 << 36 | s0 >> 28;
	return result;
}

static inline void kad_xoroshiro128plus_jump(kad_rng_t *r)
{
	static const uint64_t JUMP[] = { 0xbeac0467eba5facbULL, 0xd86b048b86aa9922ULL };
	uint64_t s0 = 0, s1 = 0;
	int i, b;
	for (i = 0; i < 2; ++i)
		for (b = 0; b < 64; b++) {
			if (JUMP[i] & 1ULL << b)
				s0 ^= r->s[0], s1 ^= r->s[1];
			kad_xoroshiro128plus_next(r);
		}
	r->s[0] = s0, r->s[1] = s1;
}

void kad_srand(void *d, uint64_t seed)
{
	kad_rng_t *r = d? (kad_rng_t*)d : &kad_rng_dat;
	r->n_gset = 0.0, r->n_iset = 0;
	r->s[0] = kad_splitmix64(seed);
	r->s[1] = kad_splitmix64(r->s[0]);
}

void *kad_rng(void)
{
	kad_rng_t *r;
	r = (kad_rng_t*)calloc(1, sizeof(kad_rng_t));
	kad_xoroshiro128plus_jump(&kad_rng_dat);
	r->s[0] = kad_rng_dat.s[0], r->s[1] = kad_rng_dat.s[1];
	return r;
}

uint64_t kad_rand(void *d) { return kad_xoroshiro128plus_next(d? (kad_rng_t*)d : &kad_rng_dat); }

double kad_drand(void *d)
{
	union { uint64_t i; double d; } u;
	u.i = 0x3FFULL << 52 | kad_xoroshiro128plus_next(d? (kad_rng_t*)d : &kad_rng_dat) >> 12;
	return u.d - 1.0;
}

// Gibbs sampling

struct gibbs_aux {
	struct { uint32_t c1:30, x:1, obs:1; } p[2];
};

static void hk_nei_gibbs1(struct hk_nei *n, struct gibbs_aux *a, float pseudo_cnt)
{
	int32_t i;
	for (i = 0; i < n->n_pairs; ++i) {
		int64_t off = n->offcnt[i] >> 16;
		int32_t cnt = n->offcnt[i] & 0xffff, j;
		struct gibbs_aux *q = &a[i];
		double s, c[2];
		c[0] = c[1] = 0.5 * pseudo_cnt, s = pseudo_cnt;
		for (j = 0; j < cnt; ++j) {
			struct hk_nei1 *n1 = &n->nei[off + j];
			c[0] += a[n1->i].p[0].x? n1->_.w : 0.0;
			c[1] += a[n1->i].p[1].x? n1->_.w : 0.0;
			s += n1->_.w;
		}
		s = 1.0 / s;
		c[0] *= s, c[1] *= s;
		if (!q->p[0].obs) q->p[0].x = (kad_drand(0) < c[0]);
		if (!q->p[1].obs) q->p[1].x = (kad_drand(0) < c[1]);
		q->p[0].c1 += q->p[0].x;
		q->p[1].c1 += q->p[1].x;
	}
}

void hk_nei_gibbs(struct hk_nei *n, struct hk_pair *pairs, int n_burnin, int n_iter, float pseudo_cnt)
{
	struct gibbs_aux *a;
	int32_t i, iter, tot = 0;

	// initialize
	a = CALLOC(struct gibbs_aux, n->n_pairs);
	for (i = 0; i < n->n_pairs; ++i) {
		a[i].p[0].x = (kad_drand(0) >= 0.5);
		a[i].p[1].x = (kad_drand(0) >= 0.5);
	}
	for (i = 0; i < n->n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (p->phase[0] >= 0) a[i].p[0].x = !!p->phase[0], a[i].p[0].obs = 1;
		if (p->phase[1] >= 0) a[i].p[1].x = !!p->phase[1], a[i].p[1].obs = 1;
	}

	// sampling, I think
	for (iter = 0; iter < n_burnin + n_iter; ++iter) {
		if (iter == n_burnin) {
			tot = 0;
			for (i = 0; i < n->n_pairs; ++i)
				a[i].p[0].c1 = a[i].p[1].c1 = 0;
		}
		hk_nei_gibbs1(n, a, pseudo_cnt);
		++tot;
		fprintf(stderr, "%d\n", iter);
	}

	// write back
	for (i = 0; i < n->n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		p->phase[0] = (int)(100.0f * a[i].p[0].c1 / tot + .499f);
		p->phase[1] = (int)(100.0f * a[i].p[1].c1 / tot + .499f);
	}
	free(a);
}

