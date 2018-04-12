#include <math.h>
#include <assert.h>
#include "hkpriv.h"

#include "ksort.h"
#define link_key(a) ((a).x)
KRADIX_SORT_INIT(link, struct hk_link, link_key, 64)

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

struct hk_link *hk_pair2link(const struct hk_opt *opt, int32_t n_pairs, struct hk_pair *pairs, int32_t *n_links_)
{
	int32_t i, j, m_links = 0, n_links = 0, k, off;
	struct hk_link *links = 0;

	// reset hk_pair::{n,offset}
	for (i = 0; i < n_pairs; ++i)
		pairs[i].n = 0, pairs[i].offset = -1;

	// find edges
	for (i = 1; i < n_pairs; ++i) {
		struct hk_pair *q = &pairs[i];
		for (j = i - 1; j >= 0; --j) {
			struct hk_pair *p = &pairs[j];
			float w;
			if (p->chr != q->chr || hk_ppos1(q) - hk_ppos1(p) > opt->max_radius)
				break;
			w = hk_pair_weight(p, q, opt->max_radius, opt->alpha, opt->beta);
			if (w > 0.0f) {
				struct hk_link *t;
				if (n_links == m_links)
					EXPAND(links, m_links);
				t = &links[n_links++];
				t->x = (uint64_t)j<<32 | i, t->w = w;
				if (n_links == m_links)
					EXPAND(links, m_links);
				t = &links[n_links++];
				t->x = (uint64_t)i<<32 | j, t->w = w;
			}
		}
	}
	fprintf(stderr, "#links: %ld\n", (long)n_links);

	// set hk_pair::{offset,n} for pairs linked to others
	rs_insertsort_link(links, links + n_links);
	fprintf(stderr, "sorted\n");
	for (off = 0, k = 1; k <= n_links; ++k) {
		if (k == n_links || links[k].x>>32 != links[k-1].x>>32) {
			struct hk_pair *p = &pairs[links[off].x>>32];
			p->offset = off;
			p->n = k - off;
			off = k;
		}
	}

	// set hk_pair::{offset,n} for orphan pairs
	int32_t n_orphans = 0;
	for (off = 0, i = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (p->offset < 0)
			p->offset = off, ++n_orphans;
		else off = p->offset;
	}
	fprintf(stderr, "#orphans/#pairs: %d/%d\n", n_orphans, n_pairs);

	*n_links_ = n_links;
	return links;
}
