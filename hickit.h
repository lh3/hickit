#ifndef HICKIT_H
#define HICKIT_H

#include <stdint.h>
#include <stdio.h>

#define HK_SUB_DELIM '!'

#ifdef __cplusplus
extern "C" {
#endif

struct hk_opt {
	int min_dist, max_seg, min_mapq;
	int min_tad_size;
	float area_weight;
	int max_radius, min_pre_link_dist;
	float alpha, beta;
};

struct hk_sdict {     // sequence dictionary
	int32_t n, m;
	char **name;
	int32_t *len;
	void *h;
};

struct hk_seg {       // a segment
	int32_t frag_id, chr, st, en;
	int8_t strand, phase;
	int16_t mapq;
};

struct hk_pair {      // a contact pair
	uint64_t chr;     // chr1<<32 | chr2
	uint64_t pos;     // pos1<<32 | pos2
	int8_t strand[2]; // strand
	int8_t phase[2];  // phase
	int32_t n;
	int64_t offset;
};

struct hk_link {
	uint64_t x;
	float w;
} __attribute__ ((__packed__));

struct hk_map {
	struct hk_sdict *d;
	int32_t n_frags, n_segs, n_pairs, n_links;
	struct hk_seg  *segs;
	struct hk_pair *pairs;
	struct hk_link *links;
};

extern int hk_verbose;

void hk_opt_init(struct hk_opt *opt);

struct hk_map *hk_map_read(const char *fn);
void hk_map_destroy(struct hk_map *m);

struct hk_pair *hk_seg2pair(int32_t n_segs, const struct hk_seg *segs, int min_dist, int max_seg, int min_mapq, int32_t *n_pairs_);
int32_t hk_pair_dedup(int n_pairs, struct hk_pair *pairs, int min_dist);
int32_t hk_pair_filter(int n_pairs, struct hk_pair *pairs, int min_dist);
void hk_pair_count(int32_t n_pairs, struct hk_pair *pairs);

struct hk_pair *hk_pair2tad(const struct hk_sdict *d, int32_t n_pairs, struct hk_pair *pairs, int min_tad_size, float area_weight, int32_t *n_tads_);
struct hk_pair *hk_pair2tad_slow(const struct hk_sdict *d, int32_t n_pairs, struct hk_pair *pairs, int max_radius, float area_weight, int32_t *n_tads_);
int32_t hk_mask_by_tad(int32_t n_tads, const struct hk_pair *tads, int32_t n_pairs, struct hk_pair *pairs);
struct hk_link *hk_pair2link(int32_t n_pairs, struct hk_pair *pairs, int max_radius, float alpha, float beta, int32_t *n_links_);

void hk_print_seg(FILE *fp, const struct hk_sdict *d, int32_t n_segs, const struct hk_seg *segs);
void hk_print_pair(FILE *fp, const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs);

#ifdef __cplusplus
}
#endif

#endif
