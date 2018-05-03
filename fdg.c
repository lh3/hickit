#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hkpriv.h"
#include "krng.h"
#include "ksort.h"
#include "kavl.h"

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

#define cy_cmp(a, b) ((a)->x[1] != (b)->x[1]? (a)->x[1] - (b)->x[1] : (a)->i - (b)->i)
KAVL_INIT(cy, struct avl_coor, head, cy_cmp)

#define cx_lt(a, b) ((a).x[0] < (b).x[0])
KSORT_INIT(cx, struct avl_coor, cx_lt)

void hk_fdg_opt_init(struct hk_fdg_opt *opt)
{
	opt->max_init = 100.0f;
	opt->n_iter = 1000;
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

static inline float fv3_norm(fvec3_t x)
{
	float s, t;
	s = sqrtf(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	t = 1.0f / s;
	x[0] *= t, x[1] *= t, x[2] *= t;
	return s;
}

static inline float fv3_sub_norm(const fvec3_t x, const fvec3_t y, fvec3_t z)
{
	z[0] = x[0] - y[0];
	z[1] = x[1] - y[1];
	z[2] = x[2] - y[2];
	return fv3_norm(z);
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

static inline void update_force(const fvec3_t *x, int32_t i, int32_t j, float k, float radius, fvec3_t *f)
{
	float dist, force;
	fvec3_t delta;
	dist = fv3_sub_norm(x[i], x[j], delta);
	force = k * (radius - dist);
	fv3_scale(force, delta);
	fv3_subfrom(delta, f[i]);
	fv3_addto(delta, f[j]);
}

void hk_fdg1(const struct hk_fdg_opt *opt, struct hk_bmap *m)
{
	const float k = 0.01f;
	const float rep_radius = 2.0f;
	const float att_radius = 1.0f;

	int32_t i, j, left;
	struct avl_coor *y, *root = 0;
	fvec3_t *f, *x = m->x;

	f = CALLOC(fvec3_t, m->n_beads);

	// attractive forces
	for (i = 0; i < m->d->n; ++i) { // backbone
		int32_t off = m->offcnt[i] >> 32;
		int32_t cnt = (int32_t)m->offcnt[i];
		for (j = 1; j < cnt; ++j)
			update_force(x, off + j - 1, off + j, k, att_radius, f);
	}
	for (i = 0; i < m->n_pairs; ++i) // contact
		update_force(x, m->pairs[i].bid[0], m->pairs[i].bid[1], k, att_radius, f);

	// repulsive forces
	y = CALLOC(struct avl_coor, m->n_beads);
	for (i = 0; i < m->n_beads; ++i) {
		y[i].x[0] = x[i][0], y[i].x[1] = x[i][1], y[i].x[2] = x[i][2];
		y[i].i = i;
	}
	ks_introsort(cx, m->n_beads, y);
	kavl_insert(cy, &root, &y[0], 0);
	for (i = 1, left = 0; i < m->n_beads; ++i) {
		// update _left_
		float x0 = y[i].x[0] - rep_radius;
		for (j = left; j < i; ++j) {
			if (y[left].x[0] >= x0) break;
			kavl_erase(cy, &root, &y[left]);
		}
		left = j;
		for (j = left; j < i; ++j) {
			if (y[j].x[1] - y[i].x[1] > rep_radius || y[j].x[1] - y[i].x[1] < -rep_radius) continue;
			if (y[j].x[2] - y[i].x[2] > rep_radius || y[j].x[2] - y[i].x[2] < -rep_radius) continue;
			update_force(x, i, j, k, rep_radius, f);
		}
		kavl_insert(cy, &root, &y[i], 0);
	}

	// free
	free(y);
	free(f);
}

void hk_fdg(const struct hk_fdg_opt *opt, struct hk_bmap *m, krng_t *rng)
{
	int32_t i;
	m->x = hk_fdg_init(rng, m->n_beads, opt->max_init);
	for (i = 0; i < opt->n_iter; ++i)
		hk_fdg1(opt, m);
}
