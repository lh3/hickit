#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "hickit.h"

#define HICKIT_VERSION "r55"

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
	fprintf(fp, "    -s INT     ignore fragments with >INT segments [%d]\n", opt->max_seg);
	fprintf(fp, "    -q INT     min mapping quality [%d]\n", opt->min_mapq);
	fprintf(fp, "    -d NUM     min distance [%d]\n", opt->min_dist);
	fprintf(fp, "    -D         don't perform duplicate removal\n");
	fprintf(fp, "  TAD calling:\n");
	fprintf(fp, "    -t         call TADs\n");
	fprintf(fp, "    -a FLOAT   area weight (smaller for bigger TADs) [%.2f]\n", opt->area_weight);
	fprintf(fp, "    -m INT     min TAD size [%d]\n", opt->min_tad_size);
	fprintf(fp, "    -M         ignore pairs contained in TADs\n");
	fprintf(fp, "  Phasing:\n");
	fprintf(fp, "    -p         perform phasing\n");
	fprintf(fp, "    -r NUM     max radius [10m]\n");
	fprintf(fp, "    -n INT     max neighbors [%d]\n", opt->max_nei);
	fprintf(fp, "    -i INT     number of iterations [%d]\n", opt->n_iter);
	fprintf(fp, "    -G         use Gibbs sampling instead of EM-like\n");
	fprintf(fp, "    -b INT     number of burn-in iterations (Gibbs only) [%d]\n", opt->n_burnin);
	fprintf(fp, "    -R INT     random seed (Gibbs only) []\n");
}

int main(int argc, char *argv[])
{
	struct hk_opt opt;
	struct hk_map *m = 0;
	int c, ret = 0, is_seg_out = 0, is_phase = 0, is_dedup = 1, is_tad_out = 0, is_gibbs = 0, sel_phased = 0, mask_tad = 0;
	int png_width = 800;
	char *fn_png = 0;

	hk_opt_init(&opt);
	while ((c = getopt(argc, argv, "o:R:SptDMGPr:v:d:s:a:m:n:fr:b:i:w:V")) >= 0) {
		if (c == 'S') is_seg_out = 1;
		else if (c == 's') opt.max_seg = atoi(optarg);
		else if (c == 'a') opt.area_weight = atof(optarg);
		else if (c == 'r') opt.max_radius = hk_parse_num(optarg);
		else if (c == 'd') opt.min_dist = hk_parse_num(optarg);
		else if (c == 'm') opt.min_tad_size = atoi(optarg);
		else if (c == 'n') opt.max_nei = atoi(optarg);
		else if (c == 'b') opt.n_burnin = atoi(optarg);
		else if (c == 'i') opt.n_iter = atoi(optarg);
		else if (c == 'f') opt.flag |= HK_OUT_PHASE;
		else if (c == 'R') kad_srand(0, atoi(optarg));
		else if (c == 'G') is_gibbs = 1;
		else if (c == 'M') mask_tad = 1;
		else if (c == 't') is_tad_out = 1;
		else if (c == 'p') is_phase = 1;
		else if (c == 'D') is_dedup = 0;
		else if (c == 'P') sel_phased = 1;
		else if (c == 'v') hk_verbose = atoi(optarg);
		else if (c == 'o') fn_png = optarg;
		else if (c == 'w') png_width = atoi(optarg);
		else if (c == 'V') {
			puts(HICKIT_VERSION);
			return 0;
		}
	}
	if (argc - optind == 0) {
		print_usage(stderr, &opt);
		return 1;
	}
	if (mask_tad && is_tad_out && hk_verbose >= 2)
		fprintf(stderr, "[W::%s] option -M is ignored\n", __func__);

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
	hk_pair_count(m->n_pairs, m->pairs);
	if (is_tad_out || mask_tad) {
		int32_t n_tads;
		struct hk_pair *tads;
		tads = hk_pair2tad(m->d, m->n_pairs, m->pairs, opt.min_tad_size, opt.area_weight, &n_tads);
		if (is_tad_out)
			hk_print_pair(stdout, opt.flag, m->d, n_tads, tads);
		else if (mask_tad)
			m->n_pairs = hk_mask_by_tad(n_tads, tads, m->n_pairs, m->pairs);
		free(tads);
		if (is_tad_out)
			goto main_return;
	}

	if (is_phase) { // phasing
		struct hk_nei *n;
		if (sel_phased)
			m->n_pairs = hk_pair_select_phased(m->n_pairs, m->pairs);
		n = hk_pair2nei(m->n_pairs, m->pairs, opt.max_radius, opt.max_nei);
		hk_nei_weight(n, opt.max_radius, opt.beta);
		if (!is_gibbs) hk_nei_phase2(n, m->pairs, opt.n_iter, opt.pseudo_cnt);
		else hk_nei_gibbs(n, m->pairs, opt.n_burnin, opt.n_iter, opt.pseudo_cnt);
		hk_nei_destroy(n);
		hk_print_pair(stdout, HK_OUT_PHASE | HK_OUT_PHASE_REAL, m->d, m->n_pairs, m->pairs);
	} else {
		if (fn_png) hk_pair_image(m->d, m->n_pairs, m->pairs, png_width, 0.1f, fn_png);
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
