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
	uint64_t pos[2];
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

struct hk_gopt {
	int min_dist, max_seg, min_mapq, max_radius;
	float alpha, beta;
};

struct hk_map *hk_map_read(const char *fn);
void hk_map_destroy(struct hk_map *m);

void hk_gopt_init(struct hk_gopt *c);
struct hk_graph *hk_graph_gen(const struct hk_map *m, const struct hk_gopt *c);
void hk_graph_destroy(struct hk_graph *g);

void hk_map_print(FILE *fp, const struct hk_map *m);
void hk_map_print_pairs(FILE *fp, const struct hk_map *m, int min_dist, int max_seg, int min_mapq);

#ifdef __cplusplus
}
#endif

#endif
