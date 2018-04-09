#include <math.h>
#include "hkpriv.h"

void hk_gopt_init(struct hk_gopt *c)
{
	c->min_dist = 500;
	c->max_seg = 3;
	c->min_mapq = 20;
	c->max_radius = 10000000;
	c->alpha = 3.0f;
	c->beta = 3.0f;
}

float hk_pair_weight(const struct hk_pair *a, const struct hk_pair *b, int32_t max, float alpha, float beta)
{
	int32_t y, z;
	float yf, zf, d;
	if (a->pos[0]>>32 != b->pos[0]>>32 || a->pos[1]>>32 != b->pos[1]>>32)
		return 0.0f;
	y = b->pos[0] > a->pos[0]? b->pos[0] - a->pos[0] : a->pos[0] - b->pos[0];
	z = b->pos[1] > a->pos[1]? b->pos[1] - a->pos[1] : a->pos[1] - b->pos[1];
	if (y > max || z > max)
		return 0.0f;
	yf = y / (float)max;
	zf = z / (float)max;
	d = 0.5f * (yf * yf + zf * zf + (alpha + alpha - 1.0f) * (yf - zf) * (yf - zf));
	return d >= 1.0f? 0.0f : expf(-beta * d);
}

struct hk_graph *hk_graph_gen(const struct hk_map *m, const struct hk_gopt *c)
{
	struct hk_graph *g;
	int32_t i, j, m_links = 0;
	g = CALLOC(struct hk_graph, 1);
	g->pairs = hk_map2pairs(m, &g->n_pairs, c->min_dist, c->max_seg, c->min_mapq);
	for (i = 1; i < g->n_pairs; ++i) {
		struct hk_pair *q = &g->pairs[i];
		for (j = i - 1; j >= 0; --j) {
			struct hk_pair *p = &g->pairs[j];
			float w;
			if (p->pos[0] - q->pos[0] > c->max_radius)
				break;
			w = hk_pair_weight(p, q, c->max_radius, c->alpha, c->beta);
			if (w > 0.0f) {
				struct hk_link *t;
				if (g->n_links == m_links)
					EXPAND(g->links, m_links);
				t = &g->links[g->n_links++];
				t->p[0] = j, t->p[1] = i, t->w = w;
				if (g->n_links == m_links)
					EXPAND(g->links, m_links);
				t = &g->links[g->n_links++];
				t->p[1] = j, t->p[0] = i, t->w = w;
			}
		}
	}
	fprintf(stderr, "%d\n", g->n_links);
	return g;
}
