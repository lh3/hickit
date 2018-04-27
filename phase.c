#include <math.h>
#include <string.h>
#include "hkpriv.h"
#include "krng.h"

static inline float dist2weight(int32_t d, int32_t max)
{
	return 1.0f / ((float)d / max + 1e-3f);
}

static float hk_pseudo_weight(int32_t max_radius)
{
	const int n = 1000;
	int i, step = max_radius / n, d;
	double sum;
	for (i = 0, d = 0, sum = 0.0; i < n; ++i, d += step) // naive integral on [0,max_radius)
		sum += dist2weight(d, max_radius);
	sum /= n;
	return sum;
}

void hk_nei_weight(struct hk_nei *n, int32_t max_radius)
{
	int32_t i;
	float coef;
	coef = 1.0 / hk_pseudo_weight(max_radius);
	for (i = 0; i < n->n_pairs; ++i) {
		int64_t off = n->offcnt[i] >> 16;
		int32_t cnt = n->offcnt[i] & 0xffff, j;
		for (j = 0; j < cnt; ++j) {
			struct hk_nei1 *n1 = &n->nei[off + j];
			n1->_.w = coef * dist2weight(n1->_.d, max_radius);
		}
	}
}

struct phase_aux {
	float p[4];
};

static inline void spacial_adj(const struct hk_pair *p1, float p[4])
{
	const float f = 1e-8f;
	int32_t d;
	float r, q[4];
	if ((int32_t)(p1->chr>>32) != (int32_t)p1->chr) return;
	if (p[0] + p[3] + 1.0f == 1.0f || p[1] + p[2] + 1.0f == 1.0f) return;
	d = hk_ppos2(p1) - hk_ppos1(p1);
	r = 1.0f - f * d;
	if (r < 0.9f || r < p[0] + p[3]) return;
	if (r > 0.99f) r = 0.99f;
	q[0] = r * p[0] / (p[0] + p[3]);
	q[3] = r * p[3] / (p[0] + p[3]);
	q[1] = (1.0f - r) * p[1] / (p[1] + p[2]);
	q[2] = (1.0f - r) * p[2] / (p[1] + p[2]);
	memcpy(p, q, 4 * sizeof(float));
}

void hk_nei_impute(struct hk_nei *n, struct hk_pair *pairs, int n_iter, float pseudo_cnt, int use_spacial)
{
	struct phase_aux *a[2], *cur, *pre, *tmp;
	int32_t i, iter;

	a[0] = CALLOC(struct phase_aux, n->n_pairs * 2);
	a[1] = a[0] + n->n_pairs;
	cur = a[0], pre = a[1];

	// initialize
	for (i = 0; i < n->n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		struct phase_aux *q = &cur[i];
		if (p->phase[0] >= 0 && p->phase[1] >= 0) { // both phase known
			q->p[p->phase[0]<<1|p->phase[1]] = 1.0f;
		} else if (p->phase[0] < 0 && p->phase[1] < 0) { // both phase unknown
			q->p[0] = q->p[1] = q->p[2] = q->p[3] = 0.25f;
		} else if (p->phase[0] >= 0) { // leg0 phase known
			q->p[p->phase[0]<<1|0] = q->p[p->phase[0]<<1|1] = 0.5f;
		} else { // p->phase[1] >= 0; leg1 phase known
			q->p[0<<1|p->phase[1]] = q->p[1<<1|p->phase[1]] = 0.5f;
		}
		if (use_spacial) spacial_adj(p, q->p);
	}
	tmp = cur, cur = pre, pre = tmp;

	// EM, I think
	for (iter = 0; iter < n_iter; ++iter) {
		for (i = 0; i < n->n_pairs; ++i) {
			struct hk_pair *p = &pairs[i];
			struct phase_aux *q = &cur[i];
			int64_t off = n->offcnt[i] >> 16;
			int32_t cnt = n->offcnt[i] & 0xffff, j;
			double c[4], s;
			c[0] = c[1] = c[2] = c[3] = 0.25 * pseudo_cnt, s = pseudo_cnt;
			for (j = 0; j < cnt; ++j) {
				struct hk_nei1 *n1 = &n->nei[off + j];
				c[0] += n1->_.w * pre[n1->i].p[0];
				c[1] += n1->_.w * pre[n1->i].p[1];
				c[2] += n1->_.w * pre[n1->i].p[2];
				c[3] += n1->_.w * pre[n1->i].p[3];
				s += n1->_.w;
			}
			for (j = 0, s = 1.0 / s; j < 4; ++j) c[j] *= s;
			if (p->phase[0] >= 0 && p->phase[1] >= 0) { // both phase known
				q->p[p->phase[0]<<1|p->phase[1]] = 1.0f;
			} else if (p->phase[0] < 0 && p->phase[1] < 0) { // both phase unknown
				q->p[0] = c[0], q->p[1] = c[1], q->p[2] = c[2], q->p[3] = c[3];
			} else if (p->phase[0] >= 0) { // leg0 phase known
				s = 1.0 / (c[p->phase[0]<<1|0] + c[p->phase[0]<<1|1]);
				q->p[p->phase[0]<<1|0] = c[p->phase[0]<<1|0] * s;
				q->p[p->phase[0]<<1|1] = c[p->phase[0]<<1|1] * s;
			} else { // leg1 phase known
				s = 1.0 / (c[0<<1|p->phase[1]] + c[1<<1|p->phase[1]]);
				q->p[0<<1|p->phase[1]] = c[0<<1|p->phase[1]] * s;
				q->p[1<<1|p->phase[1]] = c[1<<1|p->phase[1]] * s;
			}
			if (use_spacial) spacial_adj(p, q->p);
		}
		tmp = cur, cur = pre, pre = tmp;
		if (iter && iter%10 == 0 && hk_verbose >= 3)
			fprintf(stderr, "[M::%s] %d iterations done\n", __func__, iter+1);
	}

	// write back
	for (i = 0; i < n->n_pairs; ++i)
		memcpy(pairs[i]._.p4, pre[i].p, 4 * sizeof(float));
	free(a[0]);
}

struct gibbs_aux {
	struct { uint32_t c1:31, x:1; } p[4];
};

static void hk_nei_gibbs1(krng_t *r, struct hk_nei *n, struct hk_pair *pairs, struct gibbs_aux *a, float pseudo_cnt)
{
	int32_t i;
	for (i = 0; i < n->n_pairs; ++i) {
		int64_t off = n->offcnt[i] >> 16;
		int32_t cnt = n->offcnt[i] & 0xffff, which, j;
		struct gibbs_aux *q = &a[i];
		struct hk_pair *p = &pairs[i];
		double s, t, c[5];
		c[0] = c[1] = c[2] = c[3] = 0.25 * pseudo_cnt, s = pseudo_cnt;
		for (j = 0; j < cnt; ++j) {
			struct hk_nei1 *n1 = &n->nei[off + j];
			c[0] += a[n1->i].p[0].x? n1->_.w : 0.0;
			c[1] += a[n1->i].p[1].x? n1->_.w : 0.0;
			c[2] += a[n1->i].p[2].x? n1->_.w : 0.0;
			c[3] += a[n1->i].p[3].x? n1->_.w : 0.0;
			s += n1->_.w; // TODO: we can precompute this if we really care about efficiency
		}
		if (p->phase[0] >= 0 && p->phase[1] >= 0) {
			which = p->phase[0]<<1|p->phase[1];
		} else if (p->phase[0] < 0 && p->phase[1] < 0) {
			for (j = 0, s = 1.0 / s; j < 4; ++j) c[j] *= s;
			s = kr_drand_r(r);
			for (j = 0, t = 0.0; j < 3; ++j) {
				t += c[j];
				if (s < t) break;
			}
			which = j;
		} else if (p->phase[0] >= 0) {
			s = kr_drand_r(r) * (c[p->phase[0]<<1|0] + c[p->phase[0]<<1|1]);
			which = p->phase[0] << 1 | (s >= c[p->phase[0]<<1|0]);
		} else {
			s = kr_drand_r(r) * (c[0<<1|p->phase[1]] + c[1<<1|p->phase[1]]);
			which = (s >= c[0<<1|p->phase[1]]) << 1 | p->phase[1];
		}
		q->p[0].x = q->p[1].x = q->p[2].x = q->p[3].x = 0;
		q->p[which].x = 1;
		++q->p[which].c1;
	}
}

void hk_nei_gibbs(krng_t *r, struct hk_nei *n, struct hk_pair *pairs, int n_burnin, int n_iter, float pseudo_cnt)
{
	struct gibbs_aux *a;
	int32_t i, j, iter, tot;

	a = CALLOC(struct gibbs_aux, n->n_pairs);
	for (i = 0; i < n->n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		int which;
		if (p->phase[0] >= 0 && p->phase[1] >= 0) {
			which = p->phase[0] << 1 | p->phase[1];
		} else if (p->phase[0] < 0 && p->phase[1] < 0) {
			which = (int)(kr_drand_r(r) * 4.0);
		} else if (p->phase[0] >= 0) {
			which = p->phase[0] << 1 | (int)(kr_drand_r(r) * 2.0);
		} else {
			which = (int)(kr_drand_r(r) * 2.0) << 1 | p->phase[1];
		}
		a[i].p[which].x = 1;
	}
	for (iter = 0, tot = 0; iter < n_burnin + n_iter; ++iter) {
		if (iter == n_burnin)
			for (tot = 0, i = 0; i < n->n_pairs; ++i)
				a[i].p[0].c1 = a[i].p[1].c1 = a[i].p[2].c1 = a[i].p[3].c1 = 0;
		hk_nei_gibbs1(r, n, pairs, a, pseudo_cnt);
		++tot;
		fprintf(stderr, "%d\n", iter);
	}
	for (i = 0; i < n->n_pairs; ++i)
		for (j = 0; j < 4; ++j)
			pairs[i]._.p4[j] = (float)a[i].p[j].c1 / tot;
	free(a);
}

void hk_validate_holdback(krng_t *r, float ratio, int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t i;
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (p->phase[0] >= 0 || p->phase[1] >= 0) {
			if (kr_drand_r(r) < ratio) {
				if (p->phase[0] >= 0) p->phase[0] = -10 + p->phase[0];
				if (p->phase[1] >= 0) p->phase[1] = -10 + p->phase[0];
			}
		}
	}
}

void hk_validate_revert(int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t i;
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (p->phase[0] < -1 || p->phase[1] < -1) {
			if (p->phase[0] == -10 || p->phase[0] == -9) p->phase[0] += 10;
			if (p->phase[1] == -10 || p->phase[1] == -9) p->phase[1] += 10;
		}
	}
}

void hk_validate_roc(int32_t n_pairs, struct hk_pair *pairs)
{
	int32_t i, j, all[2][100], cnt[2][2][100], tot[2], sum_all[2], sum_cnt[2][2];
	memset(all, 0, sizeof(int32_t) * 100 * 2);
	memset(cnt, 0, sizeof(int32_t) * 100 * 4);
	for (i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		float max = -1.0f;
		int32_t w, max_j = -1, off_diag;
		if (p->chr>>32 != (int32_t)p->chr) off_diag = 1;
		else off_diag = !p->tad_masked;
		for (j = 0; j < 4; ++j)
			if (max < p->_.p4[j])
				max = p->_.p4[j], max_j = j;
		w = (int)(max * 100.0f);
		if (w == 100) w = 99;
		++all[off_diag][w];
		if (p->phase[0] >= -1 && p->phase[1] >= -1) continue;
		if (p->phase[0] == -10 || p->phase[0] == -9) {
			int32_t phase = p->phase[0] + 10;
			int32_t u = (max_j>>1 == phase);
			++cnt[off_diag][u][w];
		}
		if (p->phase[1] == -10 || p->phase[1] == -9) {
			int32_t phase = p->phase[1] + 10;
			int32_t u = ((max_j&1) == phase);
			++cnt[off_diag][u][w];
		}
	}
	for (i = 0, tot[0] = tot[1] = 0; i < 100; ++i)
		tot[0] += all[0][i], tot[1] += all[1][i];
	for (j = 0; j < 2; ++j)
		sum_all[j] = sum_cnt[j][0] = sum_cnt[j][1] = 0;
	for (i = 99; i > 25; --i) {
		for (j = 0; j < 2; ++j)
			sum_all[j] += all[j][i], sum_cnt[j][0] += cnt[j][0][i], sum_cnt[j][1] += cnt[j][1][i];
		printf("%.2f\t%.4f\t%.4f\t%.4f\t%.4f\n", .01f * i,
			   (float)sum_all[0] / tot[0], (sum_cnt[0][1] + .05) / (sum_cnt[0][0] + sum_cnt[0][1] + 0.1),
			   (float)sum_all[1] / tot[1], (sum_cnt[1][1] + .05) / (sum_cnt[1][0] + sum_cnt[1][1] + 0.1));
	}
	hk_validate_revert(n_pairs, pairs);
}
