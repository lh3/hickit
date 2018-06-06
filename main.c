#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "hickit.h"

#define HICKIT_VERSION "r228"

static struct option long_options_pair[] = {
	{ "out-phase",      no_argument,       0, 0 }, // 0
	{ "out-seg",        no_argument,       0, 0 }, // 1: output the segment format; mostly for testing
	{ "seed",           required_argument, 0, 'S' },
	{ "no-dedup",       no_argument,       0, 'D' },
	{ "verbose",        required_argument, 0, 0 },
	{ "select-phased",  no_argument,       0, 0 }, // 5: unused!!!
	{ "no-spacial",     no_argument,       0, 'u' },
	{ "tad-flag",       no_argument,       0, 0 }, // 7
	{ "radius",         required_argument, 0, 0 }, // 8
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

static void print_usage_pair(FILE *fp, const struct hk_popt *opt)
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
	fprintf(fp, "    -P INT        ploidy to infer sex chr phases for male (1 or 2) [2]\n");
	fprintf(fp, "  TAD calling:\n");
	fprintf(fp, "    -t            call and output TADs\n");
	fprintf(fp, "    -a FLOAT      area weight (smaller for bigger TADs) [%.2f]\n", opt->area_weight);
	fprintf(fp, "    -c INT        min TAD size [%d]\n", opt->min_tad_size);
	fprintf(fp, "  Imputation:\n");
	fprintf(fp, "    -p            impute phases with EM\n");
	fprintf(fp, "    -n INT        max neighbors within max radius [%d]\n", opt->max_nei);
	fprintf(fp, "    -k INT        number of iterations [%d]\n", opt->n_iter);
	fprintf(fp, "    -v FLOAT      fraction of phased legs held out for validation [0]\n");
	fprintf(fp, "    -u            disable the spacial heuristic (EM-only so far)\n");
	fprintf(fp, "  Miscellaneous:\n");
	fprintf(fp, "    -S INT        random seed for Gibbs sampling and validation [1]\n");
	fprintf(fp, "    --out-phase   output two 'phase' columns in .pairs\n");
}

int main_pair(int argc, char *argv[])
{
	struct hk_popt opt;
	struct hk_map *m = 0;
	int c, long_idx, ret = 0, is_seg_out = 0, is_impute = 0, is_dedup = 1, is_tad_out = 0, mask_tad = 0, use_spacial = 1;
	int cnt_radius = 0, ploidy = 2, seed = 1;
	float val_frac = -1.0f;
	krng_t rng;

	hk_popt_init(&opt);
	while ((c = getopt_long(argc, argv, "s:q:d:DP:m:ta:c:pr:n:k:v:uS:", long_options_pair, &long_idx)) >= 0) {
		if (c == 's') opt.max_seg = atoi(optarg);
		else if (c == 'q') opt.min_mapq = atoi(optarg);
		else if (c == 'd') opt.min_dist = hk_parse_num(optarg);
		else if (c == 'D') is_dedup = 0;
		else if (c == 'P') ploidy = atoi(optarg);
		else if (c == 'm') opt.min_flt_cnt = atoi(optarg);
		else if (c == 't') is_tad_out = 1;
		else if (c == 'a') opt.area_weight = atof(optarg);
		else if (c == 'c') opt.min_tad_size = atoi(optarg);
		else if (c == 'p') is_impute = 1;
		else if (c == 'r') opt.max_radius = hk_parse_num(optarg);
		else if (c == 'n') opt.max_nei = atoi(optarg);
		else if (c == 'k') opt.n_iter = atoi(optarg);
		else if (c == 'v') val_frac = atof(optarg), is_impute = 1;
		else if (c == 'u') use_spacial = 0;
		else if (c == 'S') seed = atoi(optarg);
		else if (c == 0 && long_idx == 0) opt.flag |= HK_OUT_PHASE;
		else if (c == 0 && long_idx == 1) is_seg_out = 1;
		else if (c == 0 && long_idx == 4) hk_verbose = atoi(optarg); // --verbose
		else if (c == 0 && long_idx == 7) mask_tad = 1; // --tad-flag
		else if (c == 0 && long_idx == 8) cnt_radius = hk_parse_num(optarg), opt.flag |= HK_OUT_CNT; // --radius
	}
	assert(ploidy >= 1 && ploidy <= 2);
	if (argc - optind == 0) {
		print_usage_pair(stderr, &opt);
		return 1;
	}
	kr_srand_r(&rng, seed);
	if (is_impute) mask_tad = 1;
	if (mask_tad && is_tad_out && hk_verbose >= 2)
		fprintf(stderr, "[W::%s] option --tad-flag is ignored\n", __func__);

	m = hk_map_read(argv[optind]);
	if (ploidy == 2) hk_map_phase_male_XY(m);

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
	else if (opt.min_dist > 0)
		m->n_pairs = hk_pair_filter_close_legs(m->n_pairs, m->pairs, opt.min_dist);
	if (is_dedup)
		m->n_pairs = hk_pair_dedup(m->n_pairs, m->pairs, opt.min_dist);
	if (opt.min_flt_cnt > 0)
		m->n_pairs = hk_pair_filter_isolated(m->n_pairs, m->pairs, opt.max_radius, opt.min_flt_cnt, 0.0f);
	if (cnt_radius > 0) {
		hk_pair_count_nei(m->n_pairs, m->pairs, cnt_radius);
		hk_print_pair(stdout, opt.flag, m->d, m->n_pairs, m->pairs);
		goto main_return;
	}
	if (is_tad_out || mask_tad) {
		int32_t n_tads;
		struct hk_pair *tads;
		hk_pair_count_contained(m->n_pairs, m->pairs);
		tads = hk_pair2tad(m->d, m->n_pairs, m->pairs, opt.min_tad_size, opt.area_weight, &n_tads);
		if (is_tad_out)
			hk_print_pair(stdout, opt.flag, m->d, n_tads, tads);
		else if (mask_tad)
			hk_mark_by_tad(n_tads, tads, m->n_pairs, m->pairs);
		free(tads);
		if (is_tad_out)
			goto main_return;
	}

	if (is_impute) { // phasing
		if (val_frac > 0.0f && val_frac < 1.0f)
			hk_validate_holdback(&rng, val_frac, m->n_pairs, m->pairs);
		if (hk_verbose >= 3)
			fprintf(stderr, "[M::%s] imputing phases from %d pairs\n", __func__, m->n_pairs);
		hk_impute(m->n_pairs, m->pairs, opt.max_radius, opt.min_radius, opt.max_nei, opt.n_iter, opt.pseudo_cnt, use_spacial);
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
	int c, bin_size = 1000000, min_cnt = 5, ploidy = 2, seed = 1, fdg = 0, flt_radius = 10000000, iso_radius = 1000000, out_pairs = 0;
	float phase_thres = 0.75f, drop_frac = 0.001f, max_dist = 0.0f;
	struct hk_map *m;
	struct hk_bmap *bm, *in = 0;
	struct hk_fdg_conf opt;
	char *fn_in = 0;
	krng_t rng;

	hk_fdg_conf_init(&opt);
	while ((c = getopt(argc, argv, "b:D:R:c:f:d:p:P:gi:k:r:e:n:S:x")) >= 0) {
		if (c == 'b') bin_size = hk_parse_num(optarg);
		else if (c == 'D') iso_radius = hk_parse_num(optarg);
		else if (c == 'R') flt_radius = hk_parse_num(optarg);
		else if (c == 'c') min_cnt = atoi(optarg);
		else if (c == 'f') drop_frac = atof(optarg);
		else if (c == 'd') max_dist = atof(optarg);
		else if (c == 'p') phase_thres = atof(optarg);
		else if (c == 'P') ploidy = atoi(optarg);
		else if (c == 'g') fdg = 1;
		else if (c == 'i') fn_in = optarg;
		else if (c == 'k') opt.k_rel_rep = atof(optarg);
		else if (c == 'r') opt.d_r = atof(optarg);
		else if (c == 'e') opt.step = atof(optarg);
		else if (c == 'n') opt.n_iter = atoi(optarg);
		else if (c == 'S') seed = atoi(optarg);
		else if (c == 'x') out_pairs = 1;
	}
	assert(ploidy >= 1 && ploidy <= 2);
	if (optind == argc) {
		fprintf(stderr, "Usage: hickit bin [options] <in.pairs>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Binning:\n");
		fprintf(stderr, "    -b NUM        bin size [1m]\n");
		fprintf(stderr, "    -R NUM        radius for filtering [10m]\n");
		fprintf(stderr, "    -c INT        min count [%d]\n", min_cnt);
		fprintf(stderr, "    -f FLOAT      drop FLOAT fraction of contacts [%g]\n", drop_frac);
		fprintf(stderr, "    -d FLOAT      distance to drop faraway pairs (effective with -i) [0]\n");
		fprintf(stderr, "    -p FLOAT      phase threshold [%g]\n", phase_thres);
		fprintf(stderr, "    -P INT        ploidy (1 or 2) [%d]\n", ploidy);
		fprintf(stderr, "  3D modeling with FDG:\n");
		fprintf(stderr, "    -g            perform FDG\n");
		fprintf(stderr, "    -i FILE       previous FDG output []\n");
		fprintf(stderr, "    -k FLOAT      relative repulsive stiffness [%g]\n", opt.k_rel_rep);
		fprintf(stderr, "    -r FLOAT      relative repulsive radius [%g]\n", opt.d_r);
		fprintf(stderr, "    -e FLOAT      step size [%g]\n", opt.step);
		fprintf(stderr, "    -n INT        max iteration [%d]\n", opt.n_iter);
		fprintf(stderr, "    -S INT        seed for initialization [%d]\n", seed);
		return 1;
	}
	kr_srand_r(&rng, seed);
	m = hk_map_read(argv[optind]);
	assert(m && m->pairs);
	if (ploidy == 2) {
		struct hk_map *tmp;
		tmp = hk_pair_split_phase(m, phase_thres);
		hk_map_destroy(m);
		m = tmp;
	}
	if (min_cnt > 0 || drop_frac > 0.0f)
		m->n_pairs = hk_pair_filter_isolated(m->n_pairs, m->pairs, flt_radius, min_cnt, drop_frac);
	if (iso_radius > 0)
		m->n_pairs = hk_pair_filter_isolated(m->n_pairs, m->pairs, iso_radius, 1, 0.0f);
	hk_pair_count_nei(m->n_pairs, m->pairs, flt_radius);
	if (out_pairs) {
		hk_print_pair(stdout, HK_OUT_PPROB, m->d, m->n_pairs, m->pairs);
		hk_map_destroy(m);
		return 0;
	}
	if (fn_in) in = hk_3dg_read(fn_in);
	if (in && max_dist > 1.01f)
		m->n_pairs = hk_pair_flt_3d(in, m->n_pairs, m->pairs, max_dist);
	bm = hk_bmap_gen(m->d, m->n_pairs, m->pairs, bin_size);
	hk_map_destroy(m);
	if (fdg) {
		hk_fdg(&opt, bm, in, &rng);
		hk_print_3dg(stdout, bm);
	} else hk_print_bmap(stdout, bm);
	hk_bmap_destroy(bm);
	if (in) hk_bmap_destroy(in);
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

#ifdef HAVE_GL
int main_view3d(int argc, char *argv[])
{
	int c, color_seed = 1;
	struct hk_v3d_opt opt;
	struct hk_bmap *m;
	char *hl = 0;

	hk_v3d_opt_init(&opt);
	hk_v3d_prep(&argc, argv);
	while ((c = getopt(argc, argv, "w:s:u:l:r:")) >= 0) {
		if (c == 'w') opt.width = atoi(optarg);
		else if (c == 's') color_seed = atoi(optarg);
		else if (c == 'l') opt.line_width = atof(optarg);
		else if (c == 'r') opt.bead_radius = atof(optarg);
		else if (c == 'u') hl = optarg;
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: hickit view3d [options] <in.3dg>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -w INT      viewer width [%d]\n", opt.width);
		fprintf(stderr, "  -s INT      seed for RNG to generate random colors [%d]\n", color_seed);
		fprintf(stderr, "  -l FLOAT    line width [%g]\n", opt.line_width);
		fprintf(stderr, "  -r FLOAT    bead radius [%g]\n", opt.bead_radius);
		fprintf(stderr, "  -u STR      comma-delimited list of chr to highlight []\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Key bindings:\n");
		fprintf(stderr, "  arrows: rotate\n");
		fprintf(stderr, "  [ or ]: rotate\n");
		fprintf(stderr, "  . or >: into screen\n");
		fprintf(stderr, "  , or <: out of screen\n");
		fprintf(stderr, "  b or B: toggle wire/ball mode\n");
		fprintf(stderr, "  z or Z: toggle white/black background\n");
		fprintf(stderr, "  f or F: toggle color for numerical/chromosome\n");
		fprintf(stderr, "  c or C: change colors\n");
		fprintf(stderr, "  q or Q: exit\n");
		return 1;
	}
	m = hk_3dg_read(argv[optind]);
	assert(m);
	hk_v3d_view(m, &opt, color_seed, hl);
	hk_bmap_destroy(m);
	return 0;
}
#endif

int main(int argc, char *argv[])
{
	int ret = 0;
	if (argc == 1) {
		fprintf(stderr, "Usage: hickit <command> [arguments]\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  pair       filtering, TAD calling and phase imputation\n");
		fprintf(stderr, "  bin        binning\n");
		fprintf(stderr, "  image2d    generate 2D contact map in PNG\n");
#ifdef HAVE_GL
		fprintf(stderr, "  view3d     view 3D structure\n");
#endif
		fprintf(stderr, "  version    print version number\n");
		return 1;
	}
	if (strcmp(argv[1], "pair") == 0) ret = main_pair(argc-1, argv+1);
	else if (strcmp(argv[1], "bin") == 0) ret = main_bin(argc-1, argv+1);
	else if (strcmp(argv[1], "image2d") == 0) ret = main_image2d(argc-1, argv+1);
#ifdef HAVE_GL
	else if (strcmp(argv[1], "view3d") == 0) ret = main_view3d(argc-1, argv+1);
#endif
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
