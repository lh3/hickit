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

void hk_phase(struct hk_nei *n, struct hk_pair *pairs, int n_iter)
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
			int64_t off = n->offcnt[i] >> 16;
			int32_t cnt = n->offcnt[i] & 0xffff, j;
			double c[2], s = 0.0;
			c[0] = c[1] = 0.0;
			for (j = 0; j < cnt; ++j) {
				struct hk_nei1 *n1 = &n->nei[off + j];
				c[0] += n1->_.w * pre[n1->i].phase[0];
				c[1] += n1->_.w * pre[n1->i].phase[1];
				s += n1->_.w;
			}
			cur[i].phase[0] = c[0] / s;
			cur[i].phase[1] = c[1] / s;
		}
		tmp = cur, cur = pre, pre = tmp;
	}

	// write back
	for (i = 0; i < n->n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		p->phase[0] = (int)(pre[i].phase[0] * 100.0f + .499f);
		p->phase[1] = (int)(pre[i].phase[1] * 100.0f + .499f);
	}
	free(a[0]);
}
