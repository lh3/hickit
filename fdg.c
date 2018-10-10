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

/*********************
 * Vector operations *
 *********************/

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

/******************
 * FDG parameters *
 ******************/

void hk_fdg_cal_c(struct hk_fdg_conf *opt)
{
	double t = (double)opt->d_c3 - opt->d_c2;
	assert(opt->d_c2 > 0.0 && opt->d_c3 > 0.0 && opt->d_c3 > opt->d_c2);
	opt->c_c1 = 3.0 * t;
	opt->c_c2 = t * t * t;
}

void hk_fdg_conf_init(struct hk_fdg_conf *opt)
{
	opt->target_radius = 10.0f;
	opt->n_iter = 1000;
	opt->step = 0.01f;
	opt->coef_moment = 0.9f;
	opt->max_f = 50.0f;

	opt->k_rel_rep = 0.1f;
	opt->sr_max = 1.0f;
	opt->d_b1 = 0.1f, opt->d_b2 = 1.1f;
	opt->d_c1 = 0.5f, opt->d_c2 = 1.5f, opt->d_c3 = 2.0f;
	hk_fdg_cal_c(opt);
}

/**********************
 * Basic FDG routines *
 **********************/

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

double hk_fdg_bond_dist(const struct hk_bmap *m)
{
	int32_t i, n;
	fvec3_t tmp;
	double sum;
	assert(m->x);
	for (i = 1, sum = 0.0, n = 0; i < m->n_beads; ++i)
		if (m->beads[i-1].chr == m->beads[i].chr)
			sum += fv3_sub_normalize(m->x[i-1], m->x[i], tmp), ++n;
	return sum / n;
}

int32_t hk_fdg_bead_size(const struct hk_bmap *m)
{
	int32_t mid_bead_gsize, i, *tmp;
	tmp = CALLOC(int32_t, m->n_beads);
	for (i = 0; i < m->n_beads; ++i)
		tmp[i] = m->beads[i].en - m->beads[i].st;
	mid_bead_gsize = ks_ksmall_int32_t(m->n_beads, tmp, (int)(m->n_beads * 0.5));
	free(tmp);
	return mid_bead_gsize;
}

double hk_fdg_copy_x(struct hk_bmap *dst, const struct hk_bmap *src, krng_t *rng)
{
	int32_t i;
	double avg_bb, dst_unit;
	assert(dst->d->n == src->d->n);
	if (dst->x) free(dst->x);
	avg_bb = hk_fdg_bond_dist(src);
	dst_unit = avg_bb / pow((double)dst->n_beads / src->n_beads, 1.0 / 3.0);
	dst->x = CALLOC(fvec3_t, dst->n_beads);
	for (i = 0; i < dst->n_beads; ++i) {
		struct hk_bead *pd = &dst->beads[i];
		int32_t j;
		j = hk_bmap_pos2bid(src, pd->chr, pd->st);
		if (j < src->n_beads - 1 && src->beads[j+1].chr == src->beads[j].chr) {
			float t = (float)(pd->st - src->beads[j].st) / (src->beads[j+1].st - src->beads[j].st);
			dst->x[i][0] = (1.0f - t) * src->x[j][0] + t * src->x[j+1][0] + .333f * dst_unit * (2.0 * kr_drand_r(rng) - 1.0);
			dst->x[i][1] = (1.0f - t) * src->x[j][1] + t * src->x[j+1][1] + .333f * dst_unit * (2.0 * kr_drand_r(rng) - 1.0);
			dst->x[i][2] = (1.0f - t) * src->x[j][2] + t * src->x[j+1][2] + .333f * dst_unit * (2.0 * kr_drand_r(rng) - 1.0);
		} else {
			dst->x[i][0] = src->x[j][0] + .333f * dst_unit * (2.0 * kr_drand_r(rng) - 1.0);
			dst->x[i][1] = src->x[j][1] + .333f * dst_unit * (2.0 * kr_drand_r(rng) - 1.0);
			dst->x[i][2] = src->x[j][2] + .333f * dst_unit * (2.0 * kr_drand_r(rng) - 1.0);
		}
	}
	return avg_bb;
}

/********************
 * FDG optimization *
 ********************/

#define FORCE_BACKBONE 1
#define FORCE_REPEL    2
#define FORCE_CONTACT  3

static inline float update_force(const struct hk_fdg_conf *conf, fvec3_t *x, int32_t i, int32_t j, float k, float unit, float d_scale, int force_type, fvec3_t *f, float *dist)
{
	float force, energy, r, t;
	fvec3_t delta;
	assert(i != j);
	*dist = fv3_sub_normalize(x[i], x[j], delta) / unit;
	r = *dist / d_scale;
	assert(r > 0.0f);
	if (force_type == FORCE_REPEL) {
		if (r >= 1.0f) return 0.0f;
		t = 1.0f - r;
		energy = k * t * t;
		force  = 2.0f * k * t;
	} else if (force_type == FORCE_BACKBONE) {
		if (r < conf->d_b1) {
			t = conf->d_b1 - r;
			energy = k * t * t;
			force  = 2.0f * k * t;
		} else if (r <= conf->d_b2) {
			force = energy = 0.0f;
		} else {
			t = r - conf->d_b2;
			energy = k * t * t;
			force  = -2.0f * k * t;
		}
	} else {
		if (r < conf->d_c1) {
			t = conf->d_c1 - r;
			energy = k * t * t;
			force  = 2.0f * k * t;
		} else if (r <= conf->d_c2) {
			force = energy = 0.0f;
		} else if (r <= conf->d_c3) {
			t = r - conf->d_c2;
			energy = k * t * t;
			force  = -2.0f * k * t;
		} else {
			t = r - conf->d_c2;
			energy = k * (conf->c_c1 * (r - conf->d_c3) + conf->c_c2 / t);
			force  = -k * (conf->c_c1 - conf->c_c2 / (t * t));
		}
	}
	fv3_scale(force, delta);
	fv3_addto(delta, f[i]);
	fv3_subfrom(delta, f[j]);
	return energy;
}

static double hk_fdg1(const struct hk_fdg_conf *opt, struct hk_bmap *m, khash_t(set64) *h, float unit, int max_nei, int mid_bead_gsize, float rel_rep_k, int iter, fvec3_t *x0)
{
	const double a_third = 1.0 / 3.0;
	int32_t i, j, n_y, left, n_bb = 0, n_rep = 0, n_con = 0;
	struct avl_coor *y, *root = 0;
	fvec3_t *f, *x = m->x;
	double sum = 0.0, e_bb = 0.0, e_con = 0.0, e_rep = 0.0, d_bb = 0.0, d_con = 0.0, d_rep = 0.0;
	float step, max_rep_radius, dist;

	step = opt->step * unit * pow(m->n_beads / 1500.0, a_third);
	f = CALLOC(fvec3_t, m->n_beads);

	// compute the radius of each bead
	for (i = 0; i < m->n_beads; ++i) {
		struct hk_bead *p = &m->beads[i];
		p->sr = 0.5 * pow((double)(p->en - p->st) / mid_bead_gsize, a_third); // min_bead_gsize: the median genomic length of beads
	}

	// apply attractive forces
	for (i = 0; i < m->d->n; ++i) { // backbone
		int32_t off = m->offcnt[i] >> 32;
		int32_t cnt = (int32_t)m->offcnt[i];
		for (j = 1; j < cnt; ++j) {
			int32_t bid = off + j;
			float d_opt = m->beads[bid-1].sr + m->beads[bid].sr;
			e_bb += update_force(opt, x, bid - 1, bid, 1.0f, unit, d_opt, FORCE_BACKBONE, f, &dist);
			d_bb += dist / d_opt;
			++n_bb;
		}
	}
	for (i = 0; i < m->n_pairs; ++i) { // contact
		const struct hk_bpair *p = &m->pairs[i];
		float k, d_scale;
		if (p->bid[0] == p->bid[1]) continue;
		k = p->max_nei >= max_nei? 1.0f : powf((double)p->max_nei / max_nei, a_third);
		d_scale = pow(p->n, -a_third);
		e_con += update_force(opt, x, p->bid[0], p->bid[1], k, unit, d_scale, FORCE_CONTACT, f, &dist);
		d_con += dist / d_scale;
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
	max_rep_radius = 2.0 * unit * opt->sr_max;
	kavl_insert(cy, &root, &y[0], 0);
	for (i = 1, left = 0; i < n_y; ++i) {
		struct avl_coor t, *q = &y[i];
		const struct avl_coor *p;
		kavl_itr_t(cy) itr;
		float sr_q, search_radius;
		// figure out the size of the current bead
		sr_q = m->beads[q->i].sr < opt->sr_max? m->beads[q->i].sr : opt->sr_max;
		if (sr_q < 1.0f) sr_q = 1.0f;
		search_radius = unit * (sr_q + opt->sr_max);
		// update _left_: discard beads out of range on the X axis
		float x0 = q->x[0] - max_rep_radius; // NB: has to use rep_radius here
		for (j = left; j < i; ++j) {
			if (y[j].x[0] >= x0) break;
			kavl_erase(cy, &root, &y[j], 0);
		}
		left = j;
		assert(kavl_size(head, root) == i - left);
		// traverse neighbors in 3D
		t.i = 0, t.x[0] = x0, t.x[1] = q->x[1] - search_radius, t.x[2] = q->x[2] - search_radius;
		kavl_itr_find(cy, root, &t, &itr);
		while ((p = kavl_at(&itr)) != 0) { // NB: for a small bead, repulsive forces are applied twice; for a large one, only once
			float dx, dy, dz, f2, max_d, sr_p;
			dy = p->x[1] - q->x[1];
			if (dy > search_radius) break; // out of range on the Y-axis
			sr_p = m->beads[p->i].sr < opt->sr_max? m->beads[p->i].sr : opt->sr_max;
			if (sr_p < 1.0f) sr_p = 1.0f;
			max_d = unit * (sr_p + sr_q);
			dx = p->x[0] - q->x[0];
			dz = p->x[2] - q->x[2];
			if (dx >= -max_d && dx <= max_d && dy <= max_d && dz >= -max_d && dz <= max_d) {
				khint_t k;
				k = kh_get(set64, h, (uint64_t)q->i << 32 | p->i);
				if (k == kh_end(h)) {
					f2 = update_force(opt, x, q->i, p->i, opt->k_rel_rep * rel_rep_k, unit, sr_p + sr_q, FORCE_REPEL, f, &dist);
					if (f2 > 0.0f) {
						e_rep += f2;
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
		for (j = 0; j < 3; ++j) {
			float t = x[i][j];
			x[i][j] += opt->coef_moment * (t - x0[i][j]) + f[i][j] * step;
			x0[i][j] = t;
		}
	}

	// free
	free(y);
	free(f);

	sum = sqrt(sum / m->n_beads);
	e_bb /= n_bb, e_con /= n_con, e_rep /= n_rep;
	d_bb /= n_bb, d_con /= n_con, d_rep /= n_rep;
	if (hk_verbose >= 3 && (iter + 1) % 10 == 0)
		fprintf(stderr, "[M::%s] iter:%d rep_coef:%.4f RMS_force:%.4f n_rep:%.4f energy:%.4f,%.4f,%.4f dist:%.4f,%.4f,%.4f\n",
				__func__, iter+1, rel_rep_k, sum, (float)n_rep / m->n_beads, e_bb, e_con, e_rep, d_bb, d_con, d_rep);
	return sum;
}

void hk_fdg(const struct hk_fdg_conf *opt, struct hk_bmap *m, const struct hk_bmap *src, krng_t *rng)
{
	const float alpha = 10.0f, turning = 0.333f;
	int32_t iter, i, j, absent, max_nei, *tmp, mid_bead_gsize;
	khash_t(set64) *h;
	fvec3_t *best_x, *x0;
	double best = 1e30;
	float unit;

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

	// figure out max_nei
	tmp = CALLOC(int32_t, m->n_pairs);
	for (i = 0; i < m->n_pairs; ++i)
		tmp[i] = m->pairs[i].max_nei;
	max_nei = ks_ksmall_int32_t(m->n_pairs, tmp, (int)(m->n_pairs * 0.5));
	free(tmp);

	mid_bead_gsize = hk_fdg_bead_size(m);
	if (hk_verbose >= 3) fprintf(stderr, "[M::%s] mid_bead_gsize:%d max_nei:%d\n", __func__, mid_bead_gsize, max_nei);

	// initialize
	if (src) {
		float src_dist;
		src_dist = hk_fdg_copy_x(m, src, rng);
		unit = src_dist / pow((double)m->n_beads / src->n_beads, 1.0/3.0);
	} else {
		m->x = hk_fdg_init(rng, m->n_beads, opt->target_radius);
		unit = fdg_optimal_dist(opt->target_radius, m->n_beads);
	}

	// FDG
	best_x = CALLOC(fvec3_t, m->n_beads);
	x0 = CALLOC(fvec3_t, m->n_beads);
	memcpy(x0, m->x, m->n_beads * sizeof(fvec3_t));
	for (iter = 0; iter < opt->n_iter; ++iter) {
		double s, rel_rep_k;
		rel_rep_k = (double)(iter + 1) / opt->n_iter;
		rel_rep_k = 1.0 / (1.0 + exp(-alpha * (rel_rep_k - turning)));
		//rel_rep_k = 1.0;
		s = hk_fdg1(opt, m, h, unit, max_nei, mid_bead_gsize, rel_rep_k, iter, x0);
		if (s < best) {
			memcpy(best_x, m->x, sizeof(fvec3_t) * m->n_beads);
			best = s;
		}
	}
	kh_destroy(set64, h);
	memcpy(m->x, best_x, sizeof(fvec3_t) * m->n_beads);
	free(x0);
	free(best_x);
}

int32_t hk_pair_flt_3d(const struct hk_bmap *m, int32_t n_pairs, struct hk_pair *pairs, float max_factor)
{
	int32_t i, n;
	fvec3_t tmp;
	double avg_bb;
	avg_bb = hk_fdg_bond_dist(m);
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
	int32_t i, j, n_d;
	fvec3_t center, tmp;
	float d, max_radius = 0.0f;
	double sum[3], sum_d = 0.0, sum_d2 = 0.0;

	sum[0] = sum[1] = sum[2] = 0.0;
	for (i = 0; i < m->n_beads; ++i)
		for (j = 0; j < 3; ++j) sum[j] += m->x[i][j];
	for (j = 0; j < 3; ++j) center[j] = sum[j] / m->n_beads;
	sum_d = hk_fdg_bond_dist(m);
	for (i = 1, n_d = 0; i < m->n_beads; ++i) {
		if (i > 0 && m->beads[i].chr == m->beads[i-1].chr) {
			d = fv3_sub_normalize(m->x[i], m->x[i-1], tmp);
			sum_d2 += (d - sum_d) * (d - sum_d);
			++n_d;
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
