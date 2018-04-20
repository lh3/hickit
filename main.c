#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "hickit.h"

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
	fprintf(fp, "  -s INT     ignore fragments with >INT segments [%d]\n", opt->max_seg);
	fprintf(fp, "  -q INT     min mapping quality [%d]\n", opt->min_mapq);
	fprintf(fp, "  -d NUM     min distance [%d]\n", opt->min_dist);
	fprintf(fp, "  -D         don't perform duplicate removal\n");
	fprintf(fp, "  -r NUM     max radius [10m]\n");
	fprintf(fp, "  -t         call TADs\n");
	fprintf(fp, "  -a FLOAT   area weight (smaller for bigger TADs) [%.2f]\n", opt->area_weight);
	fprintf(fp, "  -M         ignore pairs contained in TADs\n");
}

int main(int argc, char *argv[])
{
	struct hk_opt opt;
	struct hk_map *m = 0;
	int c, ret = 0, is_seg_out = 0, is_graph = 0, is_dedup = 1, is_tad_out = 0, mask_tad = 0;

	hk_opt_init(&opt);
	while ((c = getopt(argc, argv, "SgtDMr:v:d:s:a:m:n:")) >= 0) {
		if (c == 'S') is_seg_out = 1;
		else if (c == 's') opt.max_seg = atoi(optarg);
		else if (c == 'a') opt.area_weight = atof(optarg);
		else if (c == 'r') opt.max_radius = hk_parse_num(optarg);
		else if (c == 'd') opt.min_dist = hk_parse_num(optarg);
		else if (c == 'm') opt.min_tad_size = atoi(optarg);
		else if (c == 'n') opt.max_nei = atoi(optarg);
		else if (c == 'M') mask_tad = 1;
		else if (c == 't') is_tad_out = 1;
		else if (c == 'g') is_graph = 1;
		else if (c == 'D') is_dedup = 0;
		else if (c == 'v') hk_verbose = atoi(optarg);
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
			hk_print_pair(stdout, m->d, n_tads, tads);
		else if (mask_tad)
			m->n_pairs = hk_mask_by_tad(n_tads, tads, m->n_pairs, m->pairs);
		free(tads);
		if (is_tad_out)
			goto main_return;
	}

	if (is_graph) { // for testing only
		struct hk_nei *n;
		m->n_pairs = hk_pair_filter(m->n_pairs, m->pairs, opt.min_pre_link_dist);
		//m->links = hk_pair2link(m->n_pairs, m->pairs, opt.max_radius, opt.alpha, opt.beta, &m->n_links);
		n = hk_pair2nei(m->n_pairs, m->pairs, opt.max_radius, opt.max_nei);
		hk_nei_destroy(n);
	} else {
		hk_print_pair(stdout, m->d, m->n_pairs, m->pairs);
	}

main_return:
	hk_map_destroy(m);
	return ret;
}
