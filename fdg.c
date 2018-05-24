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
	opt->k_rep = 0.2f;
	opt->r_rep = 2.0f;
	opt->n_iter = 1000;
	opt->step = 0.01f;
	opt->max_f = 50.0f;
}

static float fdg_optimal_dist(float target_radius, int n_beads)
{
	return target_radius / pow(n_beads, 1.0 / 3.0);
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

#define FORCE_BACKBONE        1
#define FORCE_REPEL           2
#define FORCE_CONTACT_STRONG  3
#define FORCE_CONTACT_WEAK    4

static inline float update_force(fvec3_t *x, int32_t i, int32_t j, float k, float unit, float d_scale, int force_type, fvec3_t *f, float *dist_)
{
	const float att_min = 0.1f, att_max = 1.1f, con_min = 0.5f, con_max = 1.5f, con_cut = 2.0f;
	float dist, force, r;
	fvec3_t delta;
	assert(i != j);
	*dist_ = dist = fv3_sub_normalize(x[i], x[j], delta) / unit;
	r = dist / d_scale;
	if (force_type == FORCE_REPEL) {
		if (r >= 1.0f) return 0.0f;
		r = 1.0f - r;
		force = k * r * r;
	} else if (force_type == FORCE_BACKBONE) {
		if (r < att_min) {
			r = 1.0f - r / att_min;
			force = k * r * r;
		} else if (r < att_max) {
			force = 0.0f;
		} else {
			r = r / att_max - 1.0f;
			force = -k * r * r;
		}
	} else {
		if (r < con_min) {
			r = 1.0f - r / con_min;
			force = k * r * r;
		} else if (r < con_max) {
			force = 0.0f;
		} else if (r < con_cut || force_type == FORCE_CONTACT_STRONG) {
			r = r / con_max - 1.0f;
			force = -k * r * r;
		} else {
			r = r / con_cut - 1.0f;
			force = -k * r - k * (con_cut / con_max - 1) * (con_cut / con_max - 1);
		}
	}
	fv3_scale(force, delta);
	fv3_addto(delta, f[i]);
	fv3_subfrom(delta, f[j]);
	return fabsf(force);
}

static double hk_fdg1(const struct hk_fdg_opt *opt, struct hk_bmap *m, khash_t(set64) *h, int max_nei, int mid_dist, float rel_rep_k, int iter)
{
	const double a_third = 1.0 / 3.0;
	int32_t i, j, n_y, left, n_rep = 0, n_con = 0;
	struct avl_coor *y, *root = 0;
	fvec3_t *f, *x = m->x;
	double sum = 0.0, f_bb = 0.0, f_con = 0.0, f_rep = 0.0, d_bb = 0.0, d_con = 0.0, d_rep = 0.0;
	float step, unit, rep_radius, dist;

	unit = fdg_optimal_dist(opt->target_radius, m->n_beads);
	step = opt->step * unit;
	f = CALLOC(fvec3_t, m->n_beads);

	// apply attractive forces
	for (i = 0; i < m->d->n; ++i) { // backbone
		int32_t off = m->offcnt[i] >> 32;
		int32_t cnt = (int32_t)m->offcnt[i];
		for (j = 1; j < cnt; ++j) {
			int32_t bid = off + j;
			int32_t d = ((m->beads[bid-1].en - m->beads[bid-1].st) + (m->beads[bid].en - m->beads[bid].st)) / 2;
			float d_opt;
			d_opt = pow((double)d / mid_dist, a_third);
			f_bb += update_force(x, bid - 1, bid, 1.0f, unit, d_opt, FORCE_BACKBONE, f, &dist);
			d_bb += dist;
		}
	}
	for (i = 0; i < m->n_pairs; ++i) { // contact
		const struct hk_bpair *p = &m->pairs[i];
		float k, d_scale;
		if (p->bid[0] == p->bid[1]) continue;
		k = p->max_nei >= max_nei? 1.0f : powf((double)p->max_nei / max_nei, a_third);
		d_scale = pow(p->n, -a_third);
		f_con += update_force(x, p->bid[0], p->bid[1], k, unit, d_scale, k == 1.0f? FORCE_CONTACT_STRONG : FORCE_CONTACT_WEAK, f, &dist);
		d_con += dist;
		++n_con;
	}

	// repulsive forces: generate y[]
	y = CALLOC(struct avl_coor, m->n_beads);
	for (i = n_y = 0; i < m->n_beads; ++i) {
		struct avl_coor *p = &y[n_y++];
		p->i = i, p->x[0] = x[i][0], p->x[1] = x[i][1], p->x[2] = x[i][2];
	}
	ks_introsort(cx, n_y, y);

	// apply repulsive forces
	rep_radius = unit * opt->r_rep;
	kavl_insert(cy, &root, &y[0], 0);
	for (i = 1, left = 0; i < n_y; ++i) {
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
		// traverse neighbors in 3D
		t.i = 0, t.x[0] = x0, t.x[1] = q->x[1] - rep_radius, t.x[2] = q->x[2] - rep_radius;
		kavl_itr_find(cy, root, &t, &itr);
		while ((p = kavl_at(&itr)) != 0) {
			float dz, f2;
			if (p->x[1] - q->x[1] > rep_radius) break; // out of range on the Y-axis
			dz = p->x[2] - q->x[2];
			if (dz >= -rep_radius && dz <= rep_radius) {
				khint_t k;
				k = kh_get(set64, h, (uint64_t)q->i << 32 | p->i);
				if (k == kh_end(h)) {
					f2 = update_force(x, q->i, p->i, opt->k_rep * rel_rep_k, unit, opt->r_rep, FORCE_REPEL, f, &dist);
					if (f2 > 0.0f) {
						f_rep += f2;
						d_rep += dist;
						++n_rep;
					}
				}
			}
			if (!kavl_itr_next(cy, &itr)) break;
		}
		kavl_insert(cy, &root, q, 0);
	}

	// update coordinate
	for (i = 0; i < m->n_beads; ++i) {
		float t;
		assert(!isnan(sum));
		t = fv3_L2(f[i]);
		sum += t; // TODO: check if precision is good enough
		t = sqrtf(t);
		if (t > opt->max_f)
			for (j = 0; j < 3; ++j)
				f[i][j] = f[i][j] / t * opt->max_f;
		fv3_axpy(step, f[i], m->x[i]);
	}

	// free
	free(y);
	free(f);

	sum = sqrt(sum / m->n_beads);
	f_bb /= m->n_beads - m->d->n, f_con /= n_con, f_rep /= n_rep;
	d_bb /= m->n_beads - m->d->n, d_con /= n_con, d_rep /= n_rep;
	if (hk_verbose >= 3 && (iter + 1) % 10 == 0)
		fprintf(stderr, "[M::%s] iter:%d rep_coef:%.4f RMS_force:%.4f n_rep:%.4f forces:%.4f,%.4f,%.4f dist:%.4f,%.4f,%.4f\n",
				__func__, iter+1, rel_rep_k, sum, (float)n_rep / m->n_beads, f_bb, f_con, f_rep, d_bb, d_con, d_rep);
	return sum;
}

void hk_fdg(const struct hk_fdg_opt *opt, struct hk_bmap *m, krng_t *rng)
{
	const float alpha = 10.0f;
	int32_t iter, i, j, absent, max_nei, *tmp, mid_dist;
	khash_t(set64) *h;
	fvec3_t *best_x;
	double best = 1e30;

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

	// median backbone size
	tmp = CALLOC(int32_t, m->n_beads);
	for (i = 0; i < m->n_beads; ++i)
		tmp[i] = m->beads[i].en - m->beads[i].st;
	mid_dist = ks_ksmall_int32_t(m->n_beads, tmp, (int)(m->n_beads * 0.5));
	free(tmp);

	// figure out max_nei
	tmp = CALLOC(int32_t, m->n_pairs);
	for (i = 0; i < m->n_pairs; ++i)
		tmp[i] = m->pairs[i].max_nei;
	max_nei = ks_ksmall_int32_t(m->n_pairs, tmp, (int)(m->n_pairs * 0.5));
	free(tmp);
	if (hk_verbose >= 3) fprintf(stderr, "[M::%s] mid_dist:%d max_nei:%d\n", __func__, mid_dist, max_nei);

	// FDG
	if (m->x == 0) m->x = hk_fdg_init(rng, m->n_beads, opt->target_radius);
	best_x = CALLOC(fvec3_t, m->n_beads);
	for (iter = 0; iter < opt->n_iter; ++iter) {
		double s, rel_rep_k;
		rel_rep_k = 0.5 + atan(alpha * (2.0 * (iter + 1) / opt->n_iter - 1)) / M_PI;
		rel_rep_k = 1.0;
		s = hk_fdg1(opt, m, h, max_nei, mid_dist, rel_rep_k, iter);
		if (s < best) {
			memcpy(best_x, m->x, sizeof(fvec3_t) * m->n_beads);
			best = s;
		}
	}
	kh_destroy(set64, h);
	memcpy(m->x, best_x, sizeof(fvec3_t) * m->n_beads);
	free(best_x);
}

int32_t hk_pair_flt_3d(const struct hk_bmap *m, int32_t n_pairs, struct hk_pair *pairs, float max_factor)
{
	int32_t i, n;
	fvec3_t tmp;
	double sum, avg_bb;

	assert(m->x);
	for (i = 1, sum = 0.0, n = 0; i < m->n_beads; ++i)
		if (m->beads[i-1].chr == m->beads[i].chr)
			sum += fv3_sub_normalize(m->x[i-1], m->x[i], tmp), ++n;
	avg_bb = sum / n;

	for (i = n = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		int bid[2];
		float d;
		bid[0] = hk_bmap_pos2bid(m, p->chr>>32,      hk_ppos1(p));
		bid[1] = hk_bmap_pos2bid(m, (int32_t)p->chr, hk_ppos2(p));
		d = fv3_sub_normalize(m->x[bid[0]], m->x[bid[1]], tmp);
		if (d < avg_bb * max_factor) pairs[n++] = pairs[i];
		// else fprintf(stderr, "%s\t%d\t%s\t%d\t%f\t%d\n", m->d->name[p->chr>>32], hk_ppos1(p), m->d->name[(int32_t)p->chr], hk_ppos2(p), d / avg_bb, p->n);
	}
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] filtered %d (out of %d) pairs\n", __func__, n_pairs - n, n_pairs);
	return n;
}

void hk_fdg_normalize(struct hk_bmap *m)
{
	int32_t i, j, n_d = 0;
	fvec3_t center, tmp;
	float d, max_radius = 0.0f;
	double sum[3], sum_d = 0.0, sum_d2 = 0.0;

	sum[0] = sum[1] = sum[2] = 0.0;
	for (i = 0; i < m->n_beads; ++i) {
		for (j = 0; j < 3; ++j) sum[j] += m->x[i][j];
		if (i > 0 && m->beads[i].chr == m->beads[i-1].chr) {
			d = fv3_sub_normalize(m->x[i], m->x[i-1], tmp);
			sum_d += d, ++n_d;
		}
	}
	for (j = 0; j < 3; ++j) center[j] = sum[j] / m->n_beads;
	sum_d /= n_d;
	for (i = 1; i < m->n_beads; ++i) {
		if (i > 0 && m->beads[i].chr == m->beads[i-1].chr) {
			d = fv3_sub_normalize(m->x[i], m->x[i-1], tmp);
			sum_d2 += (d - sum_d) * (d - sum_d);
		}
	}
	sum_d2 = sqrt(sum_d2 / n_d) / sum_d;
	for (i = 0; i < m->n_beads; ++i) {
		d = fv3_sub_normalize(m->x[i], center, tmp);
		max_radius = max_radius > d? max_radius : d;
	}
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] avg = %f; cv = %.4f; max_radius = %.4f\n", __func__, sum_d, sum_d2, max_radius);
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] center: (%.4f,%.4f,%.4f)\n", __func__, center[0], center[1], center[2]);
	for (i = 0; i < m->n_beads; ++i)
		for (j = 0; j < 3; ++j)
			m->x[i][j] = (m->x[i][j] - center[j]) / max_radius;
}
