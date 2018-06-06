#include <getopt.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "hickit.h"

static struct option long_options[] = {
	{ "min-leg-dist",   required_argument, 0, 0 },   // 0
	{ "max-seg",        required_argument, 0, 0 },   // 1
	{ "min-mapq",       required_argument, 0, 0 },   // 2
	{ "fn-pairs",       required_argument, 0, 'o' }, // 3
	{ "fn-seg",         required_argument, 0, 0 },   // 4
	{ "area-weight",    required_argument, 0, 0 },   // 5
	{ "tad-min-size",   required_argument, 0, 0 },   // 6
	{ "fn-tads",        required_argument, 0, 'T' }, // 7
	{ "seed",           required_argument, 0, 0 },   // 8
	{ "keep-dup",       no_argument,       0, 0 },   // 9
	{ "imput-nei",      required_argument, 0, 0 },   // 10
	{ "val-frac",       required_argument, 0, 0 },   // 11
	{ "impute",         no_argument,       0, 0 },   // 12
	{ "fn-val",         required_argument, 0, 0 },   // 13
	{ "fn-png",         required_argument, 0, 0 },   // 14
	{ "png-no-dim",     no_argument,       0, 0 },   // 15
	{ 0, 0, 0, 0}
};

static inline int64_t hk_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(optarg, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

int main_stream(int argc, char *argv[])
{
	int c, long_idx, has_options = 0;
	FILE *fp;
	struct hk_map *m = 0;
	int n_tads = 0;
	struct hk_pair *tads = 0;
	krng_t rng;

	// general parameters
	int ploidy = 2, max_iter = 1000, radius = 10000000, width = 780;
	float phase_thres = 0.75f;
	// immediate input filters
	int max_seg = 3, min_mapq = 20, min_leg_dist = 1000, dedup = 1;
	// TAD calling parameters
	float tad_area_weight = 5.0f;
	int tad_min_size = 10;
	// imputation parameters
	int imput_max_nei = 50, imput_min_radius = 50000;
	float imput_val_frac = 0.1f, imput_pseudo_cnt = 0.4f;
	// PNG generation
	int png_no_dim = 0;

	kr_srand_r(&rng, 1); // initialize the RNG with seed 1

	while ((c = getopt_long(argc, argv, "i:o:r:m:T:P:n:w:p:", long_options, &long_idx)) >= 0) {
		has_options = 1;
		if (c == 'i') {
			if (m) hk_map_destroy(m);
			m = hk_map_read(optarg);
			assert(m);
			if (ploidy == 2) hk_map_phase_male_XY(m);
			if (m->pairs == 0 && m->segs)
				m->pairs = hk_seg2pair(m->n_segs, m->segs, min_leg_dist, max_seg, min_mapq, &m->n_pairs);
			else if (m->pairs && min_leg_dist > 0)
				m->n_pairs = hk_pair_filter_close_legs(m->n_pairs, m->pairs, min_leg_dist);
			if (dedup)
				m->n_pairs = hk_pair_dedup(m->n_pairs, m->pairs, min_leg_dist);
		} else if (c == 'o') {
			fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
			assert(fp);
			hk_print_pair(fp, m->cols, m->d, m->n_pairs, m->pairs);
			if (fp != stdout) fclose(fp);
		} else if (c == 'r') {
			radius = hk_parse_num(optarg);
			assert(radius > 0);
		} else if (c == 'n') {
			max_iter = atoi(optarg);
			assert(max_iter > 0);
		} else if (c == 'p') {
			phase_thres = atof(optarg);
			assert(phase_thres > 0.0f && phase_thres <= 1.0f);
		} else if (c == 'P') {
			ploidy = atoi(optarg);
			assert(ploidy == 1 || ploidy == 2);
		} else if (c == 'w') {
			width = atoi(optarg);
			assert(width > 0);
		} else if (c == 'T') {
			assert(m && m->pairs);
			tads = hk_pair2tad(m->d, m->n_pairs, m->pairs, tad_min_size, tad_area_weight, &n_tads);
			m->cols |= 1<<7;
			fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
			hk_print_pair(fp, 1<<7, m->d, n_tads, tads);
			if (fp != stdout) fclose(fp);
		} else if (c == 'm') {
			int min_flt_cnt;
			min_flt_cnt = atoi(optarg);
			if (min_flt_cnt > 0) {
				m->n_pairs = hk_pair_filter_isolated(m->n_pairs, m->pairs, radius, min_flt_cnt, 0.0f);
				m->cols |= 1<<6;
			}
		} else if (c == 0) {
			if (long_idx == 0) min_leg_dist = hk_parse_num(optarg); // --min-leg-dist
			else if (long_idx ==  1) max_seg = atoi(optarg); // --max-seg
			else if (long_idx ==  2) min_mapq = atoi(optarg); // --min-mapq
			else if (long_idx ==  5) tad_area_weight = atof(optarg); // --area-weight
			else if (long_idx ==  6) tad_min_size = atoi(optarg); // --tad-min-size
			else if (long_idx ==  7) kr_srand_r(&rng, atol(optarg)); // --seed
			else if (long_idx ==  9) dedup = 0; // --keep-dup
			else if (long_idx == 10) imput_max_nei = atoi(optarg); // --imput-nei
			else if (long_idx == 11) imput_val_frac = atof(optarg); // --val-frac
			else if (long_idx == 15) png_no_dim = 1; // --png-no-dim
			else if (long_idx ==  4) { // --fn-seg
				assert(m && m->segs);
				fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
				assert(fp);
				hk_print_seg(fp, m->d, m->n_segs, m->segs);
				if (fp != stdout) fclose(fp);
			} else if (long_idx == 12) { // --impute
				m->cols |= 0x3c;
				hk_impute(m->n_pairs, m->pairs, radius, imput_min_radius, imput_max_nei, max_iter, imput_pseudo_cnt, 1);
			} else if (long_idx == 13) { // --fn-val
				fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
				assert(fp);
				if (tads == 0) tads = hk_pair2tad(m->d, m->n_pairs, m->pairs, tad_min_size, tad_area_weight, &n_tads);
				m->cols |= 1<<7;
				hk_mark_by_tad(n_tads, tads, m->n_pairs, m->pairs);
				hk_validate_holdback(&rng, imput_val_frac, m->n_pairs, m->pairs);
				hk_impute(m->n_pairs, m->pairs, radius, imput_min_radius, imput_max_nei, max_iter, imput_pseudo_cnt, 1);
				hk_validate_roc(fp, m->n_pairs, m->pairs);
				if (fp != stdout) fclose(fp);
			} else if (long_idx == 14) { // --fn-png
				hk_pair_image(m->d, m->n_pairs, m->pairs, width, phase_thres, png_no_dim, optarg);
			}
		}
	}

	if (!has_options) {
		FILE *fp = stderr;
		fprintf(fp, "Usage: hickit stream [options]\n");
		fprintf(fp, "Options:\n");
		fprintf(fp, "  Actions:\n");
		fprintf(fp, "    -i FILE             read .pairs or .seg FILE and apply input filters []\n");
		fprintf(fp, "    -o FILE             output .pairs to FILE []\n");
		fprintf(fp, "    -m INT              filter out isolated pairs [0]\n");
		fprintf(fp, "    -T FILE             call TADs and write to FILE []\n");
		fprintf(fp, "    --impute            impute missing phases\n");
		fprintf(fp, "    --fn-val=FILE       save validation to FILE []\n");
		fprintf(fp, "    --fn-png=FILE       write 2D contact map to FILE in PNG []\n");
		fprintf(fp, "  General settings:\n");
		fprintf(fp, "    --seed=INT          random number seed [1]\n");
		fprintf(fp, "    -P INT              ploidy: 1 or 2 [%d]\n", ploidy);
		fprintf(fp, "    -r NUM              radius [10m]\n");
		fprintf(fp, "    -n INT              max iterations [%d]\n", max_iter);
		fprintf(fp, "    -w INT              image or viewer width [%d]\n", width);
		fprintf(fp, "    -p FLOAT            prob. threshold for a contact considered phased [%g]\n", phase_thres);
		fprintf(fp, "  Input filters:\n");
		fprintf(fp, "    --max-seg=NUM       ignore fragments with >INT segments [%d]\n", max_seg);
		fprintf(fp, "    --min-mapq=NUM      min mapping quality [%d]\n", min_mapq);
		fprintf(fp, "    --min-leg-dist=NUM  min base-pair distance between the two legs [%d]\n", min_leg_dist);
		fprintf(fp, "    --keep-dup          don't filter potential duplicates\n");
		fprintf(fp, "  TAD calling:\n");
		fprintf(fp, "    --area-weight=FLOAT area weight [%g]\n", tad_area_weight);
		fprintf(fp, "    --tad-min-size=INT  min TAD size [%d]\n", tad_min_size);
		fprintf(fp, "  Imputation:\n");
		fprintf(fp, "    --imput-nei=INT     max neighbors [%d]\n", imput_max_nei);
		fprintf(fp, "    --val-frac=FLOAT    fraction to hold out for validation [%g]\n", imput_val_frac);
		return 1;
	}

	free(tads);
	hk_map_destroy(m);

	return 0;
}
