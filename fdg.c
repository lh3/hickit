#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "hkpriv.h"
#include "krng.h"
#include "ksort.h"
#include "kavl.h"
#include "khash.h"

KHASH_INIT(set64, uint64_t, char, 0, hash64, kh_int64_hash_equal)

struct fdg_coor {
	fvec3_t x;
	int32_t i;
};

// define AVL node
struct avl_coor {
	fvec3_t x;
	int32_t i;
	KAVL_HEAD(struct avl_coor) head;
};

#define cy_cmp(a, b) ((a)->x[1] != (b)->x[1]? ((a)->x[1] > (b)->x[1]) - ((a)->x[1] < (b)->x[1]) : (a)->i - (b)->i)
KAVL_INIT(cy, struct avl_coor, head, cy_cmp)

#define cx_lt(a, b) ((a).x[0] < (b).x[0])
KSORT_INIT(cx, struct avl_coor, cx_lt)

void hk_fdg_opt_init(struct hk_fdg_opt *opt)
{
	opt->target_radius = 10.0f;
	opt->k_rep = 1.0f;
	opt->r_rep = 1.0f;
	opt->n_iter = 1000;
	opt->step = 0.01f;
}

static float fdg_optimal_dist(float target_radius, int n_beads)
{
	double v;
	v = 4.0 / 3.0 * M_PI * pow(target_radius, 3.0) / n_beads;
	return pow(v, 1.0 / 3.0);
}

fvec3_t *hk_fdg_init(krng_t *rng, int n_beads, float max)
{
	int32_t i;
	fvec3_t *x;
	x = CALLOC(fvec3_t, n_beads);
	for (i = 0; i < n_beads; ++i) {
		x[i][0] = max * (2.0 * kr_drand_r(rng) - 1.0);
		x[i][1] = max * (2.0 * kr_drand_r(rng) - 1.0);
		x[i][2] = max * (2.0 * kr_drand_r(rng) - 1.0);
	}
	return x;
}

static inline float fv3_L2(fvec3_t x)
{
	return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
}

static inline float fv3_normalize(fvec3_t x)
{
	float s, t;
	s = sqrtf(fv3_L2(x));
	t = 1.0f / s;
	x[0] *= t, x[1] *= t, x[2] *= t;
	return s;
}

static inline float fv3_sub_normalize(const fvec3_t x, const fvec3_t y, fvec3_t z)
{
	z[0] = x[0] - y[0];
	z[1] = x[1] - y[1];
	z[2] = x[2] - y[2];
	return fv3_normalize(z);
}

static inline void fv3_addto(const fvec3_t x, fvec3_t y)
{
	y[0] += x[0], y[1] += x[1], y[2] += x[2];
}

static inline void fv3_subfrom(const fvec3_t x, fvec3_t y)
{
	y[0] -= x[0], y[1] -= x[1], y[2] -= x[2];
}

static inline void fv3_scale(float a, fvec3_t x)
{
	x[0] *= a, x[1] *= a, x[2] *= a;
}

static inline void fv3_axpy(float a, const fvec3_t x, fvec3_t y)
{
	y[0] += a * x[0], y[1] += a * x[1], y[2] += a * x[2];
}

static inline void update_force(const fvec3_t *x, int32_t i, int32_t j, float k, float radius, int repel, fvec3_t *f)
{
	float dist, force;
	fvec3_t delta;
	assert(i != j);
	dist = fv3_sub_normalize(x[i], x[j], delta);
	if (repel && dist >= radius) return;
	force = k * (radius - dist);
	fv3_scale(force, delta);
	fv3_addto(delta, f[i]);
	fv3_subfrom(delta, f[j]);
}

static double hk_fdg1(const struct hk_fdg_opt *opt, struct hk_bmap *m, khash_t(set64) *h)
{
	int32_t i, j, left;
	struct avl_coor *y, *root = 0;
	fvec3_t *f, *x = m->x;
	double sum = 0.0;
	float att_radius, rep_radius;

	att_radius = fdg_optimal_dist(opt->target_radius, m->n_beads);
	rep_radius = att_radius * opt->r_rep;
	f = CALLOC(fvec3_t, m->n_beads);

	// attractive forces
	for (i = 0; i < m->d->n; ++i) { // backbone
		int32_t off = m->offcnt[i] >> 32;
		int32_t cnt = (int32_t)m->offcnt[i];
		for (j = 1; j < cnt; ++j)
			update_force(x, off + j - 1, off + j, 1.0f, att_radius, 0, f);
	}
	for (i = 0; i < m->n_pairs; ++i) { // contact
		const struct hk_bpair *p = &m->pairs[i];
		if (p->bid[0] != p->bid[1])
			update_force(x, p->bid[0], p->bid[1], 1.0f, att_radius, 0, f);
	}

	// repulsive forces
	y = CALLOC(struct avl_coor, m->n_beads);
	for (i = 0; i < m->n_beads; ++i) {
		y[i].x[0] = x[i][0], y[i].x[1] = x[i][1], y[i].x[2] = x[i][2];
		y[i].i = i;
	}
	ks_introsort(cx, m->n_beads, y);
	kavl_insert(cy, &root, &y[0], 0);
	for (i = 1, left = 0; i < m->n_beads; ++i) {
		struct avl_coor t, *q = &y[i];
		const struct avl_coor *p;
		kavl_itr_t(cy) itr;
		// update _left_
		float x0 = q->x[0] - rep_radius;
		for (j = left; j < i; ++j) {
			if (y[j].x[0] >= x0) break;
			kavl_erase(cy, &root, &y[j]);
		}
		left = j;
		assert(kavl_size(head, root) == i - left);
		t.i = 0, t.x[0] = x0, t.x[1] = q->x[1] - rep_radius, t.x[2] = q->x[2] - rep_radius;
		kavl_itr_find(cy, root, &t, &itr);
		while ((p = kavl_at(&itr)) != 0) {
			float dz;
			if (p->x[1] - q->x[1] > rep_radius) break;
			dz = p->x[2] - q->x[2];
			if (dz >= -rep_radius && dz <= rep_radius) {
				khint_t k;
				k = kh_get(set64, h, (uint64_t)q->i << 32 | p->i);
				if (k == kh_end(h))
					update_force(x, q->i, p->i, opt->k_rep, rep_radius, 1, f);
			}
			if (!kavl_itr_next(cy, &itr)) break;
		}
		kavl_insert(cy, &root, q, 0);
	}

	// update coordinate
	for (i = 0; i < m->n_beads; ++i) {
		sum += fv3_L2(f[i]); // TODO: check if precision is good enough
		fv3_axpy(opt->step, f[i], m->x[i]);
	}

	// free
	free(y);
	free(f);
	return sum / m->n_beads;
}

void hk_fdg(const struct hk_fdg_opt *opt, struct hk_bmap *m, krng_t *rng)
{
	int32_t iter, i, j, absent;
	khash_t(set64) *h;

	// collect attractive pairs
	h = kh_init(set64);
	for (i = 0; i < m->d->n; ++i) {
		int32_t off = m->offcnt[i] >> 32;
		int32_t cnt = (int32_t)m->offcnt[i];
		for (j = 1; j < cnt; ++j) {
			kh_put(set64, h, (uint64_t)(off + j - 1) << 32 | (off + j), &absent);
			kh_put(set64, h, (uint64_t)(off + j) << 32 | (off + j - 1), &absent);
		}
	}
	for (i = 0; i < m->n_pairs; ++i) { // contact
		const struct hk_bpair *p = &m->pairs[i];
		kh_put(set64, h, (uint64_t)p->bid[0] << 32 | p->bid[1], &absent);
		kh_put(set64, h, (uint64_t)p->bid[1] << 32 | p->bid[0], &absent);
	}

	// FDG
	m->x = hk_fdg_init(rng, m->n_beads, opt->target_radius * 10.0f);
	for (iter = 0; iter < opt->n_iter; ++iter) {
		double s;
		s = hk_fdg1(opt, m, h);
		if (iter && iter%10 == 0 && hk_verbose >= 3)
			fprintf(stderr, "[M::%s] %d iterations done (%.4f)\n", __func__, iter+1, s);
	}
	kh_destroy(set64, h);
}
