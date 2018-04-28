#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "hickit.h"

#define HICKIT_VERSION "r83"

static struct option long_options[] = {
	{ "version",        no_argument,       0, 0 },
	{ "out-seg",        no_argument,       0, 0 }, // 1: output the segment format; mostly for testing
	{ "seed",           required_argument, 0, 'S' },
	{ "no-dedup",       no_argument,       0, 'D' },
	{ "verbose",        required_argument, 0, 0 },
	{ "select-phased",  no_argument,       0, 0 }, // 5: only imput pairs containing at least one phased leg
	{ "no-spacial",     no_argument,       0, 'u' },
	{ "tad-flag",       no_argument,       0, 0 }, // 7
	{ "out-phase",      no_argument,       0, 0 }, // 8: output two "phase" columns in .pairs
	{ "sort",           no_argument,       0, 0 }, // 9
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

static void print_usage(FILE *fp, const struct hk_opt *opt)
{
	fprintf(fp, "Usage: hickit [options] <in.seg>|<in.pairs>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Input filtering:\n");
	fprintf(fp, "    -s INT        ignore fragments with >INT segments [%d]\n", opt->max_seg);
	fprintf(fp, "    -q INT        min mapping quality [%d]\n", opt->min_mapq);
	fprintf(fp, "    -d NUM        min distance between legs [%d]\n", opt->min_dist);
	fprintf(fp, "    -D            don't remove duplicates\n");
	fprintf(fp, "    --sort        force to re-sort pairs\n");
	fprintf(fp, "  TAD calling:\n");
	fprintf(fp, "    -t            call and output TADs\n");
	fprintf(fp, "    -a FLOAT      area weight (smaller for bigger TADs) [%.2f]\n", opt->area_weight);
	fprintf(fp, "    -c INT        min TAD size [%d]\n", opt->min_tad_size);
	fprintf(fp, "    --tad-flag    flag pairs contained in TADs (forced with -p/-G)\n");
	fprintf(fp, "  Imputation:\n");
	fprintf(fp, "    -p            impute phases with EM\n");
	fprintf(fp, "    -r NUM        max radius [10m]\n");
	fprintf(fp, "    -n INT        max neighbors [%d]\n", opt->max_nei);
	fprintf(fp, "    -k INT        number of iterations [%d]\n", opt->n_iter);
	fprintf(fp, "    -G            impute with Gibbs sampling (NOT recommended)\n");
	fprintf(fp, "    -B INT        number of burn-in iterations (Gibbs only) [%d]\n", opt->n_burnin);
	fprintf(fp, "    -v FLOAT      fraction of phased legs held out for validation [0]\n");
	fprintf(fp, "    -u            disable the spacial heuristic (EM-only so far)\n");
	fprintf(fp, "  Image generation:\n");
	fprintf(fp, "    -I FILE       write PNG to FILE (no output to stdout) []\n");
	fprintf(fp, "    -w INT        width of the image [800]\n");
	fprintf(fp, "    -P FLOAT      probability threshold for a dot considered phased [0.7]\n");
	fprintf(fp, "  Miscellaneous:\n");
	fprintf(fp, "    -S INT        random seed for Gibbs sampling and validation [1]\n");
	fprintf(fp, "    --out-phase   output two 'phase' columns in .pairs\n");
	fprintf(fp, "    --version     version number\n");
}

int main(int argc, char *argv[])
{
	struct hk_opt opt;
	struct hk_map *m = 0;
	int c, long_idx, ret = 0, is_seg_out = 0, is_impute = 0, is_dedup = 1, is_tad_out = 0, is_gibbs = 0, sel_phased = 0, mask_tad = 0, use_spacial = 1;
	int assume_sorted = 1, seed = 1, png_width = 800;
	float phase_thres = 0.7f, val_frac = -1.0f;
	char *fn_png = 0;
	krng_t rng;

	hk_opt_init(&opt);
	while ((c = getopt_long(argc, argv, "s:q:d:Dta:c:pr:n:k:GB:v:I:w:P:S:", long_options, &long_idx)) >= 0) {
		if (c == 's') opt.max_seg = atoi(optarg);
		else if (c == 'q') opt.min_mapq = atoi(optarg);
		else if (c == 'd') opt.min_dist = hk_parse_num(optarg);
		else if (c == 'D') is_dedup = 0;
		else if (c == 't') is_tad_out = 1;
		else if (c == 'a') opt.area_weight = atof(optarg);
		else if (c == 'c') opt.min_tad_size = atoi(optarg);
		else if (c == 'p') is_impute = 1;
		else if (c == 'r') opt.max_radius = hk_parse_num(optarg);
		else if (c == 'n') opt.max_nei = atoi(optarg);
		else if (c == 'k') opt.n_iter = atoi(optarg);
		else if (c == 'G') is_gibbs = is_impute = 1;
		else if (c == 'B') opt.n_burnin = atoi(optarg);
		else if (c == 'v') val_frac = atof(optarg), is_impute = 1;
		else if (c == 'u') use_spacial = 0;
		else if (c == 'I') fn_png = optarg;
		else if (c == 'w') png_width = atoi(optarg);
		else if (c == 'P') phase_thres = atof(optarg);
		else if (c == 'S') seed = atoi(optarg);
		else if (c == 0 && long_idx == 1) is_seg_out = 1;
		else if (c == 0 && long_idx == 4) hk_verbose = atoi(optarg);
		else if (c == 0 && long_idx == 5) sel_phased = 1;
		else if (c == 0 && long_idx == 7) mask_tad = 1;
		else if (c == 0 && long_idx == 8) opt.flag |= HK_OUT_PHASE;
		else if (c == 0 && long_idx == 9) assume_sorted = 0;
		else if (c == 0 && long_idx == 0) {
			puts(HICKIT_VERSION);
			return 0;
		}
	}
	if (argc - optind == 0) {
		print_usage(stderr, &opt);
		return 1;
	}
	kr_srand_r(&rng, seed);
	if (is_impute) mask_tad = 1;
	if (mask_tad && is_tad_out && hk_verbose >= 2)
		fprintf(stderr, "[W::%s] option --tad-flag is ignored\n", __func__);

	m = hk_map_read(argv[optind]);

	if (is_seg_out) {
		if (m->segs) {
			hk_print_seg(stdout, m->d, m->n_segs, m->segs);
		} else {
			if (hk_verbose >= 1)
				fprintf(stderr, "[E::%s] the input is not in the seg format\n", __func__);
			ret = 1;
		}
		goto main_return;
	}

	if (m->pairs == 0 && m->segs)
		m->pairs = hk_seg2pair(m->n_segs, m->segs, opt.min_dist, opt.max_seg, opt.min_mapq, &m->n_pairs);
	else if (m->segs == 0 && m->pairs && !assume_sorted)
		hk_pair_sort(m->n_pairs, m->pairs);
	if (is_dedup)
		m->n_pairs = hk_pair_dedup(m->n_pairs, m->pairs, opt.min_dist);
	if (is_tad_out || mask_tad) {
		int32_t n_tads;
		struct hk_pair *tads;
		hk_pair_count(m->n_pairs, m->pairs);
		tads = hk_pair2tad(m->d, m->n_pairs, m->pairs, opt.min_tad_size, opt.area_weight, &n_tads);
		if (is_tad_out)
			hk_print_pair(stdout, opt.flag, m->d, n_tads, tads);
		else if (mask_tad)
			hk_mask_by_tad(n_tads, tads, m->n_pairs, m->pairs);
		free(tads);
		if (is_tad_out)
			goto main_return;
	}

	if (is_impute) { // phasing
		struct hk_nei *n;
		if (sel_phased)
			m->n_pairs = hk_pair_select_phased(m->n_pairs, m->pairs);
		n = hk_pair2nei(m->n_pairs, m->pairs, opt.max_radius, opt.max_nei);
		hk_nei_weight(n, m->pairs, opt.min_radius, opt.max_radius);
		if (val_frac > 0.0f && val_frac < 1.0f)
			hk_validate_holdback(&rng, val_frac, m->n_pairs, m->pairs);
		if (!is_gibbs) hk_nei_impute(n, m->pairs, opt.n_iter, opt.pseudo_cnt, use_spacial);
		else hk_nei_gibbs(&rng, n, m->pairs, opt.n_burnin, opt.n_iter, opt.pseudo_cnt);
		hk_nei_destroy(n);
		if (val_frac > 0.0f && val_frac < 1.0f)
			hk_validate_roc(m->n_pairs, m->pairs);
		else
			hk_print_pair(stdout, HK_OUT_P4, m->d, m->n_pairs, m->pairs);
	} else {
		if (fn_png) hk_pair_image(m->d, m->n_pairs, m->pairs, png_width, phase_thres, fn_png);
		else hk_print_pair(stdout, opt.flag, m->d, m->n_pairs, m->pairs);
	}

main_return:
	hk_map_destroy(m);
	if (hk_verbose >= 3) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, HICKIT_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fputc('\n', stderr);
	}
	return ret;
}
