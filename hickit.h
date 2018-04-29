#ifndef HICKIT_H
#define HICKIT_H

#include <stdint.h>
#include <stdio.h>
#include "krng.h"

#define HK_SUB_DELIM '!'

#define HK_OUT_PHASE        0x1
#define HK_OUT_P4           0x2

#ifdef __cplusplus
extern "C" {
#endif

struct hk_opt {
	int flag;
	int min_dist, max_seg, min_mapq, min_flt_cnt;
	int min_tad_size;
	float area_weight;
	int min_radius, max_radius, max_nei;
	float pseudo_cnt;
	int n_burnin, n_iter;
	int n_multi_ploidy;
	float phase_thres;
	int min_bin_cnt;
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
	uint32_t n:31, tad_masked:1;
	union {
		float p4[4];
	} _;
};

struct hk_nei1 {
	union {
		int32_t d;
		float w;
	} _;
	int32_t i;
};

struct hk_nei {
	int32_t n_pairs;
	uint64_t *offcnt;
	struct hk_nei1 *nei;
};

struct hk_map {
	struct hk_sdict *d;
	int32_t n_frags, n_segs, n_pairs;
	struct hk_seg  *segs;
	struct hk_pair *pairs;
};

struct hk_bpair {
	uint64_t chr;
	uint64_t pos;
	uint64_t end;
	int32_t n;
	float p;
};

struct hk_bmap {
	int32_t ploidy, n_full;
	struct hk_sdict *d;
	int32_t n_pairs;
	struct hk_bpair *pairs;
};

extern int hk_verbose;

void hk_opt_init(struct hk_opt *opt);

struct hk_sdict *hk_sd_dup(const struct hk_sdict *d, int ploidy, int n_full);

struct hk_map *hk_map_read(const char *fn);
void hk_map_destroy(struct hk_map *m);

struct hk_pair *hk_seg2pair(int32_t n_segs, const struct hk_seg *segs, int min_dist, int max_seg, int min_mapq, int32_t *n_pairs_);
int hk_pair_is_sorted(int32_t n_pairs, const struct hk_pair *pairs);
void hk_pair_sort(int32_t n_pairs, struct hk_pair *pairs);
int32_t hk_pair_dedup(int n_pairs, struct hk_pair *pairs, int min_dist);
int32_t hk_pair_filter(int32_t n_pairs, struct hk_pair *pairs, int32_t max_radius, int32_t min_cnt);
int32_t hk_pair_select_phased(int n_pairs, struct hk_pair *pairs);
void hk_pair_count(int32_t n_pairs, struct hk_pair *pairs);

struct hk_pair *hk_pair2tad(const struct hk_sdict *d, int32_t n_pairs, struct hk_pair *pairs, int min_tad_size, float area_weight, int32_t *n_tads_);
void hk_mask_by_tad(int32_t n_tads, const struct hk_pair *tads, int32_t n_pairs, struct hk_pair *pairs);

struct hk_nei *hk_pair2nei(int n_pairs, const struct hk_pair *pairs, int max_radius, int max_nei);
void hk_nei_weight(struct hk_nei *n, const struct hk_pair *pairs, int min_radius, int32_t max_radius);
void hk_nei_impute(struct hk_nei *n, struct hk_pair *pairs, int n_iter, float pseudo_cnt, int use_spacial);
void hk_nei_gibbs(krng_t *r, struct hk_nei *n, struct hk_pair *pairs, int n_burnin, int n_iter, float pseudo_cnt);
void hk_nei_destroy(struct hk_nei *n);

void hk_validate_holdback(krng_t *r, float ratio, int32_t n_pairs, struct hk_pair *pairs);
void hk_validate_roc(int32_t n_pairs, struct hk_pair *pairs);

struct hk_bmap *hk_bmap_gen2(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int size, int n_full, int min_cnt, float phase_thres);
void hk_bmap_destroy(struct hk_bmap *m);

void hk_print_seg(FILE *fp, const struct hk_sdict *d, int32_t n_segs, const struct hk_seg *segs);
void hk_print_pair(FILE *fp, int flag, const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs);

void hk_pair_image(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int w, float phase_thres, int no_grad, const char *fn);

#ifdef __cplusplus
}
#endif

#endif
