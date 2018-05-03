#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "hickit.h"

#define HICKIT_VERSION "r97"

static struct option long_options_pair[] = {
	{ "out-phase",      no_argument,       0, 0 }, // 0
	{ "out-seg",        no_argument,       0, 0 }, // 1: output the segment format; mostly for testing
	{ "seed",           required_argument, 0, 'S' },
	{ "no-dedup",       no_argument,       0, 'D' },
	{ "verbose",        required_argument, 0, 0 },
	{ "select-phased",  no_argument,       0, 0 }, // 5: only imput pairs containing at least one phased leg
	{ "no-spacial",     no_argument,       0, 'u' },
	{ "tad-flag",       no_argument,       0, 0 }, // 7
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

static void print_usage_pair(FILE *fp, const struct hk_opt *opt)
{
	fprintf(fp, "Usage: hickit pair [options] <in.seg>|<in.pairs>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Input filtering:\n");
	fprintf(fp, "    -s INT        ignore fragments with >INT segments [%d]\n", opt->max_seg);
	fprintf(fp, "    -q INT        min mapping quality [%d]\n", opt->min_mapq);
	fprintf(fp, "    -d NUM        min distance between legs [%d]\n", opt->min_dist);
	fprintf(fp, "    -r NUM        max radius (affecting imputation as well) [10m]\n");
	fprintf(fp, "    -m INT        min count within max radius [%d]\n", opt->min_flt_cnt);
	fprintf(fp, "    -D            don't remove duplicates\n");
	fprintf(fp, "  TAD calling:\n");
	fprintf(fp, "    -t            call and output TADs\n");
	fprintf(fp, "    -a FLOAT      area weight (smaller for bigger TADs) [%.2f]\n", opt->area_weight);
	fprintf(fp, "    -c INT        min TAD size [%d]\n", opt->min_tad_size);
	fprintf(fp, "  Imputation:\n");
	fprintf(fp, "    -p            impute phases with EM\n");
	fprintf(fp, "    -n INT        max neighbors within max radius [%d]\n", opt->max_nei);
	fprintf(fp, "    -k INT        number of iterations [%d]\n", opt->n_iter);
	fprintf(fp, "    -G            impute with Gibbs sampling (NOT recommended)\n");
	fprintf(fp, "    -B INT        number of burn-in iterations (Gibbs only) [%d]\n", opt->n_burnin);
	fprintf(fp, "    -v FLOAT      fraction of phased legs held out for validation [0]\n");
	fprintf(fp, "    -u            disable the spacial heuristic (EM-only so far)\n");
	fprintf(fp, "  Miscellaneous:\n");
	fprintf(fp, "    -S INT        random seed for Gibbs sampling and validation [1]\n");
	fprintf(fp, "    --out-phase   output two 'phase' columns in .pairs\n");
}

int main_pair(int argc, char *argv[])
{
	struct hk_opt opt;
	struct hk_map *m = 0;
	int c, long_idx, ret = 0, is_seg_out = 0, is_impute = 0, is_dedup = 1, is_tad_out = 0, is_gibbs = 0, sel_phased = 0, mask_tad = 0, use_spacial = 1;
	int seed = 1;
	float val_frac = -1.0f;
	krng_t rng;

	hk_opt_init(&opt);
	while ((c = getopt_long(argc, argv, "s:q:d:Dm:ta:c:pr:n:k:GB:v:uS:", long_options_pair, &long_idx)) >= 0) {
		if (c == 's') opt.max_seg = atoi(optarg);
		else if (c == 'q') opt.min_mapq = atoi(optarg);
		else if (c == 'd') opt.min_dist = hk_parse_num(optarg);
		else if (c == 'D') is_dedup = 0;
		else if (c == 'm') opt.min_flt_cnt = atoi(optarg);
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
		else if (c == 'S') seed = atoi(optarg);
		else if (c == 0 && long_idx == 0) opt.flag |= HK_OUT_PHASE;
		else if (c == 0 && long_idx == 1) is_seg_out = 1;
		else if (c == 0 && long_idx == 4) hk_verbose = atoi(optarg);
		else if (c == 0 && long_idx == 5) sel_phased = 1;
		else if (c == 0 && long_idx == 7) mask_tad = 1;
	}
	if (argc - optind == 0) {
		print_usage_pair(stderr, &opt);
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
	if (is_dedup)
		m->n_pairs = hk_pair_dedup(m->n_pairs, m->pairs, opt.min_dist);
	if (opt.min_flt_cnt > 0)
		m->n_pairs = hk_pair_filter(m->n_pairs, m->pairs, opt.max_radius, opt.min_flt_cnt);
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
		hk_print_pair(stdout, opt.flag, m->d, m->n_pairs, m->pairs);
	}

main_return:
	hk_map_destroy(m);
	return ret;
}

int main_bin(int argc, char *argv[])
{
	int c, bin_size = 1000000, min_cnt = 1, ploidy = 2, n_multi_ploidy = 23, seed = 1, fdg = 0;
	float phase_thres = 0.7f;
	struct hk_map *m;
	struct hk_bmap *bm;
	struct hk_fdg_opt opt;
	krng_t rng;

	hk_fdg_opt_init(&opt);
	while ((c = getopt(argc, argv, "c:b:p:P:f:gk:r:e:n:s:")) >= 0) {
		if (c == 'c') min_cnt = atoi(optarg);
		else if (c == 'b') bin_size = hk_parse_num(optarg);
		else if (c == 'p') phase_thres = atof(optarg);
		else if (c == 'P') ploidy = atoi(optarg);
		else if (c == 'f') n_multi_ploidy = atoi(optarg);
		else if (c == 'g') fdg = 1;
		else if (c == 'k') opt.k_rep = atof(optarg);
		else if (c == 'r') opt.r_rep = atof(optarg);
		else if (c == 'e') opt.step = atof(optarg);
		else if (c == 'n') opt.n_iter = atoi(optarg);
		else if (c == 's') seed = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: hickit bin [options] <in.pairs>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Binning:\n");
		fprintf(stderr, "    -b NUM        bin size [1m]\n");
		fprintf(stderr, "    -c INT        min count [%d]\n", min_cnt);
		fprintf(stderr, "    -p FLOAT      phase threshold [%g]\n", phase_thres);
		fprintf(stderr, "    -P INT        ploidy [%d]\n", ploidy);
		fprintf(stderr, "    -f INT        first INT chr have multiple ploidy [%d]\n", n_multi_ploidy);
		fprintf(stderr, "  FDG:\n");
		fprintf(stderr, "    -g            perform FDG\n");
		fprintf(stderr, "    -k FLOAT      relative repulsive stiffness [%g]\n", opt.k_rep);
		fprintf(stderr, "    -r FLOAT      relative repulsive radius [%g]\n", opt.r_rep);
		fprintf(stderr, "    -e FLOAT      step size [%g]\n", opt.step);
		fprintf(stderr, "    -n INT        max iteration [%d]\n", opt.n_iter);
		fprintf(stderr, "    -s INT        seed for initialization [%d]\n", seed);
		return 1;
	}
	kr_srand_r(&rng, seed);
	m = hk_map_read(argv[optind]);
	assert(m && m->pairs);
	bm = hk_bmap_gen(m->d, m->n_pairs, m->pairs, bin_size, min_cnt);
	if (ploidy > 1 && 1) {
		struct hk_bmap *bm2;
		bm2 = hk_bmap_dup(bm, ploidy, n_multi_ploidy, min_cnt, phase_thres);
		hk_bmap_destroy(bm);
		bm = bm2;
	}
	if (fdg) {
		hk_fdg(&opt, bm, &rng);
		hk_print_fdg(stdout, bm);
	} else hk_print_bmap(stdout, bm);
	hk_bmap_destroy(bm);
	hk_map_destroy(m);
	return 0;
}

int main_image2d(int argc, char *argv[])
{
	int c, width = 800, no_dim = 0;
	float phase_thres = 0.7f;
	char *fn_png = 0;
	struct hk_map *m;
	while ((c = getopt(argc, argv, "w:o:p:d")) >= 0) {
		if (c == 'w') width = atoi(optarg);
		else if (c == 'o') fn_png = optarg;
		else if (c == 'p') phase_thres = atof(optarg);
		else if (c == 'd') no_dim = 1;
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: hickit image2d [options] <in.pairs>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -o FILE      PNG file name (required) []\n");
		fprintf(stderr, "  -w INT       image width [%d]\n", width);
		fprintf(stderr, "  -p FLOAT     probability threshold for a dot considered phased [%g]\n", phase_thres);
		fprintf(stderr, "  -d           don't dim by counts\n");
		return 1;
	}
	if (fn_png == 0) {
		fprintf(stderr, "[E::%s] option -o is required (to avoid corrupting the terminal)\n", __func__);
		return 1;
	}
	m = hk_map_read(argv[optind]);
	assert(m && m->pairs);
	hk_pair_image(m->d, m->n_pairs, m->pairs, width, phase_thres, no_dim, fn_png);
	hk_map_destroy(m);
	return 0;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	if (argc == 1) {
		fprintf(stderr, "Usage: hickit <command> [arguments]\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  pair       filtering, TAD calling and phase imputation\n");
		fprintf(stderr, "  bin        binning\n");
		fprintf(stderr, "  image2d    generate 2D contact map in PNG\n");
		fprintf(stderr, "  version    print version number\n");
		return 1;
	}
	if (strcmp(argv[1], "pair") == 0) ret = main_pair(argc-1, argv+1);
	else if (strcmp(argv[1], "bin") == 0) ret = main_bin(argc-1, argv+1);
	else if (strcmp(argv[1], "image2d") == 0) ret = main_image2d(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", HICKIT_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
		ret = 1;
	}
	if (hk_verbose >= 3 && ret == 0) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, HICKIT_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fputc('\n', stderr);
	}
	return ret;
}
