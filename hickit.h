#ifndef HICKIT_H
#define HICKIT_H

#include <stdint.h>
#include <stdio.h>

#define HK_SUB_DELIM '!'

#ifdef __cplusplus
extern "C" {
#endif

struct hk_sdict {
	int32_t n, m;
	char **name;
	int32_t *len;
	void *h;
};

struct hk_seg {
	int32_t frag_id, chr, st, en;
	int8_t strand, phase;
	int16_t mapq;
};

struct hk_map {
	struct hk_sdict *d;
	int32_t n_frag, m_seg, n_seg;
	struct hk_seg *seg;
};

struct hk_pair {
	uint64_t chr, pos;
	int8_t phase[2];
	int8_t rel_strand;
	int32_t n_nei;
	int64_t offset;
};

struct hk_link {
	uint64_t x;
	float w;
} __attribute__ ((__packed__));

struct hk_graph {
	int32_t n_pairs;
	int64_t n_links;
	struct hk_pair *pairs;
	struct hk_link *links;
};

struct hk_opt {
	int min_dist, max_seg, min_mapq, max_radius;
	float alpha, beta;
};

struct hk_map *hk_map_read(const char *fn);
void hk_map_destroy(struct hk_map *m);

struct hk_pair *hk_map2pairs(const struct hk_map *m, int32_t *_n_pairs, int min_dist, int max_seg, int min_mapq);
int32_t hk_pair_dedup(int n_pairs, struct hk_pair *pairs, int min_dist);
void hk_pair_print(FILE *fp, const struct hk_map *m, int32_t n_pairs, const struct hk_pair *pairs);

void hk_opt_init(struct hk_opt *opt);
struct hk_graph *hk_graph_gen(const struct hk_map *m, const struct hk_opt *opt);
void hk_graph_destroy(struct hk_graph *g);

void hk_map_print(FILE *fp, const struct hk_map *m);

#ifdef __cplusplus
}
#endif

#endif
