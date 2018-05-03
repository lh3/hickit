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
	int32_t bid[2];
	int32_t n;
	int8_t phase[2];
	float p;
};

struct hk_bead {
	int32_t chr, st, en;
};

typedef float fvec3_t[3];

struct hk_bmap {
	int32_t ploidy, n_full, n_beads, n_pairs;
	struct hk_sdict *d;
	struct hk_bead *beads;
	uint64_t *offcnt; // index into beads
	struct hk_bpair *pairs;
	fvec3_t *x;
};

struct hk_fdg_opt {
	float max_init;
	int n_iter;
	float k;
	float rep_radius, att_radius;
};

extern int hk_verbose;

void hk_opt_init(struct hk_opt *opt);

struct hk_sdict *hk_sd_dup(const struct hk_sdict *d, int ploidy, int n_full);

struct hk_map *hk_map_read(const char *fn);
void hk_map_destroy(struct hk_map *m);

struct hk_pair *hk_seg2pair(int32_t n_segs, const struct hk_seg *segs, int min_dist, int max_seg, int min_mapq, int32_t *n_pairs_);
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

struct hk_bmap *hk_bmap_gen(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int size, int min_cnt);
struct hk_bmap *hk_bmap_dup(const struct hk_bmap *m0, int ploidy, int n_with_homo, int min_cnt, float phase_thres);
void hk_bmap_destroy(struct hk_bmap *m);

void hk_fdg_opt_init(struct hk_fdg_opt *opt);
void hk_fdg(const struct hk_fdg_opt *opt, struct hk_bmap *m, krng_t *rng);

void hk_print_seg(FILE *fp, const struct hk_sdict *d, int32_t n_segs, const struct hk_seg *segs);
void hk_print_pair(FILE *fp, int flag, const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs);
void hk_print_bmap(FILE *fp, const struct hk_bmap *m);
void hk_print_fdg(FILE *fp, const struct hk_bmap *m);

void hk_pair_image(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int w, float phase_thres, int no_grad, const char *fn);

#ifdef __cplusplus
}
#endif

#endif
