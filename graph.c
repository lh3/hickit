#include <math.h>
#include <assert.h>
#include "hkpriv.h"

#include "ksort.h"
#define link_key(a) ((a).x)
KRADIX_SORT_INIT(link, struct hk_link, link_key, 64)

#define ppos1(p) ((int32_t)((p)->pos>>32))
#define ppos2(p) ((int32_t)(p)->pos)

void hk_gopt_init(struct hk_gopt *c)
{
	c->min_dist = 500;
	c->max_seg = 3;
	c->min_mapq = 20;
	c->max_radius = 10000000;
	c->alpha = 3.0f;
	c->beta = 3.0f;
}

int32_t hk_pair_dedup(int n_pairs, struct hk_pair *pairs, int min_dist)
{
	int32_t i, j, n;
	for (i = n = 1; i < n_pairs; ++i) {
		struct hk_pair *q = &pairs[i];
		int32_t to_skip = 0;
		for (j = i - 1; j >= 0; --j) {
			struct hk_pair *p = &pairs[j];
			int32_t d;
			if (p->chr != q->chr || ppos1(q) - ppos1(p) >= min_dist)
				break;
			d = ppos2(q) - ppos2(p);
			if (d < min_dist && d > -min_dist) {
				to_skip = 1;
				break;
			}
		}
		if (!to_skip) pairs[n++] = pairs[i];
	}
	return n;
}

float hk_pair_weight(const struct hk_pair *a, const struct hk_pair *b, int32_t max, float alpha, float beta)
{
	int32_t y, z;
	float yf, zf, d;
	if (a->chr != b->chr) return 0.0f;
	y = ppos1(b) > ppos1(a)? ppos1(b) - ppos1(a) : ppos1(a) - ppos1(b);
	z = ppos2(b) > ppos2(a)? ppos2(b) - ppos2(a) : ppos2(a) - ppos2(b);
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
	int32_t i, j;
	int64_t m_links = 0, k, off;

	g = CALLOC(struct hk_graph, 1);
	g->pairs = hk_map2pairs(m, &g->n_pairs, c->min_dist, c->max_seg, c->min_mapq);
	g->n_pairs = hk_pair_dedup(g->n_pairs, g->pairs, c->min_dist);

	// find edges
	for (i = 1; i < g->n_pairs; ++i) {
		struct hk_pair *q = &g->pairs[i];
		for (j = i - 1; j >= 0; --j) {
			struct hk_pair *p = &g->pairs[j];
			float w;
			if (p->chr != q->chr || ppos1(q) - ppos1(p) > c->max_radius)
				break;
			w = hk_pair_weight(p, q, c->max_radius, c->alpha, c->beta);
			if (w > 0.0f) {
				struct hk_link *t;
				if (g->n_links == m_links)
					EXPAND(g->links, m_links);
				t = &g->links[g->n_links++];
				t->x = (uint64_t)j<<32 | i, t->w = w;
				if (g->n_links == m_links)
					EXPAND(g->links, m_links);
				t = &g->links[g->n_links++];
				t->x = (uint64_t)i<<32 | j, t->w = w;
			}
		}
	}
	fprintf(stderr, "#links: %ld\n", (long)g->n_links);

	// set hk_pair::{offset,n_nei} for pairs linked to others
	rs_insertsort_link(g->links, g->links + g->n_links);
	for (off = 0, k = 1; k <= g->n_links; ++k) {
		if (k == g->n_links || g->links[k].x>>32 != g->links[k-1].x>>32) {
			struct hk_pair *p = &g->pairs[g->links[off].x>>32];
			p->offset = off;
			p->n_nei = k - off;
			off = k;
		}
	}

	// set hk_pair::{offset,n_nei} for orphan pairs
	int32_t n_orphans = 0;
	for (off = 0, i = 0; i < g->n_pairs; ++i) {
		struct hk_pair *p = &g->pairs[i];
		if (p->offset < 0)
			p->offset = off, ++n_orphans;
		else off = p->offset;
	}
	fprintf(stderr, "#orphans/#pairs: %d/%d\n", n_orphans, g->n_pairs);
	return g;
}

void hk_graph_destroy(struct hk_graph *g)
{
	free(g->links);
	free(g->pairs);
	free(g);
}
