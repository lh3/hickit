#include "hkpriv.h"

/*********
 * Pairs *
 *********/

#include "ksort.h"
#define pair_lt(a, b) ((a).chr < (b).chr || ((a).chr == (b).chr && (a).pos < (b).pos))
KSORT_INIT(pair, struct hk_pair, pair_lt)

struct hk_pair *hk_map2pairs(const struct hk_map *m, int32_t *_n_pairs, int min_dist, int max_seg, int min_mapq)
{
	int32_t m_pairs = 0, n_pairs = 0, i, j, st;
	struct hk_pair *pairs = 0;
	for (st = 0, i = 1; i <= m->n_seg; ++i) {
		if (i == m->n_seg || m->seg[i].frag_id != m->seg[i-1].frag_id) {
			if (i - st <= max_seg) {
				for (j = st + 1; j < i; ++j) {
					struct hk_pair *p;
					struct hk_seg *s = &m->seg[j], *t = &m->seg[j-1];
					if (s->mapq < min_mapq || t->mapq < min_mapq)
						continue; // mapping quality too low
					if (n_pairs == m_pairs)
						EXPAND(pairs, m_pairs);
					p = &pairs[n_pairs++];
					if (t->chr < s->chr || (t->chr == s->chr && t->en < s->st)) {
						p->chr = (uint64_t)t->chr << 32 | s->chr;
						p->pos = (uint64_t)t->en  << 32 | s->st;
						p->phase[0] = t->phase, p->phase[1] = s->phase;
					} else {
						p->chr = (uint64_t)s->chr << 32 | t->chr;
						p->pos = (uint64_t)s->st  << 32 | t->en;
						p->phase[0] = s->phase, p->phase[1] = t->phase;
					}
					p->rel_strand = t->strand * s->strand;
					p->n = 0, p->offset = -1;
					if (p->rel_strand >= 0 && t->chr == s->chr && (int32_t)p->pos - (int32_t)(p->pos>>32) < min_dist) {
						--n_pairs;
						continue;
					}
				}
			}
			st = i;
		}
	}
	ks_introsort_pair(n_pairs, pairs);
	*_n_pairs = n_pairs;
	return pairs;
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
			if (p->chr != q->chr || hk_ppos1(q) - hk_ppos1(p) >= min_dist)
				break;
			d = hk_ppos2(q) - hk_ppos2(p);
			if (d < min_dist && d > -min_dist) {
				to_skip = 1;
				break;
			}
		}
		if (!to_skip) pairs[n++] = pairs[i];
	}
	if (hk_verbose >= 3)
		fprintf(stderr, "[M::%s] duplicate rate: %.2f%% = %d / %d\n", __func__,
				100.0 * (n_pairs - n) / n_pairs, n_pairs - n, n_pairs);
	return n;
}

void hk_pair_print(FILE *fp, const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs)
{
	int32_t i;
	fprintf(fp, "## pairs format v1.0\n");
	fprintf(fp, "#sorted: chr1-chr2-pos1-pos2\n");
	fprintf(fp, "#shape: upper triangle\n");
	hk_print_chr(fp, d);
	fprintf(fp, "#columns: readID chr1 pos1 chr2 pos2\n");
	for (i = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		fprintf(fp, ".\t%s\t%d\t%s\t%d\n", d->name[p->chr>>32], (int32_t)(p->pos>>32),
				d->name[(int32_t)p->chr], (int32_t)p->pos);
	}
}

/***************
 * TAD calling *
 ***************/

struct hk_pair *hk_tad_call1(int32_t n_pairs, struct hk_pair *pairs, int max_radius, float alpha, int32_t *_n_tads)
{
	int32_t i, n_a, min = INT32_MAX, max = INT32_MIN;
	struct hk_pair *a;
	float area_tot;

	// generate a[]
	for (i = n_a = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		if (hk_ppos2(p) - hk_ppos1(p) <= max_radius)
			++n_a;
	}
	a = MALLOC(struct hk_pair, n_a);
	for (i = n_a = 0; i < n_pairs; ++i) {
		struct hk_pair *p = &pairs[i];
		int32_t p1 = hk_ppos1(p), p2 = hk_ppos2(p);
		if (p2 - p1 <= max_radius) {
			a[n_a] = *p;
			a[n_a++].n = 0;
			min = min < p1? min : p1;
			max = max > p2? max : p2;
		}
	}
	area_tot = 1e-12 * (max - min + 0.5 * max_radius) * max_radius; // #pairs per square-million
	fprintf(stderr, "%d => %d; [%d,%d]; %f\n", n_pairs, n_a, min, max, n_a / area_tot);

	// count
	for (i = 1; i < n_a; ++i) {
		struct hk_pair *q = &a[i];
		int32_t j, q1 = hk_ppos1(q), q2 = hk_ppos2(q);
		for (j = i - 1; j >= 0; --j) {
			struct hk_pair *p = &a[j];
			int32_t p1 = hk_ppos1(p), p2 = hk_ppos2(p);
			if (q2 - p1 > max_radius) break;
			if (q1 != p1 && q2 < p2) ++p->n;
		}
	}

	// DP
	for (i = n_a - 1; i >= 0; --i) {
		struct hk_pair *p = &a[i];
		int32_t j, p1 = hk_ppos1(p), p2 = hk_ppos2(p);
		for (j = i + 1; j < n_a; ++j) {
			struct hk_pair *q = &a[j];
			int32_t q1 = hk_ppos1(q), q2;
			if (q1 - p1 > max_radius) break;
			if (p2 >= q1) continue;
			q2 = hk_ppos2(q);
		}
	}
	free(a);
	return 0;
}

struct hk_pair *hk_tad_call(const struct hk_sdict *d, int32_t n_pairs, struct hk_pair *pairs, int max_radius, float alpha, int32_t *_n_tads)
{
	int32_t i, st, n;
	struct hk_pair *tads = 0;
	for (st = 0, i = 1; i <= n_pairs; ++i) {
		if (i == n_pairs || pairs[i].chr != pairs[i-1].chr) {
			if (pairs[st].chr>>32 == (int32_t)pairs[st].chr) {
				fprintf(stderr, "===> chr: %s <===\n", d->name[(int32_t)pairs[st].chr]);
				hk_tad_call1(i - st, &pairs[st], max_radius, alpha, &n);
			}
			st = i;
		}
	}
	return tads;
}
