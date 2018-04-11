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
	fprintf(fp, "  -a FLOAT   area weight [%.2f]\n", opt->area_weight);
	fprintf(fp, "  -M         ignore pairs contained in TADs\n");
}

int main(int argc, char *argv[])
{
	struct hk_opt opt;
	struct hk_map *m = 0;
	struct hk_graph *g = 0;
	int c, is_pairs = 1, is_graph = 0, is_dedup = 1, is_tad = 0, mask_tad = 0;

	hk_opt_init(&opt);
	while ((c = getopt(argc, argv, "SgtDMr:v:d:s:a:")) >= 0) {
		if (c == 'S') is_pairs = 0;
		else if (c == 's') opt.max_seg = atoi(optarg);
		else if (c == 'a') opt.area_weight = atof(optarg);
		else if (c == 'r') opt.max_radius = hk_parse_num(optarg);
		else if (c == 'd') opt.min_dist = hk_parse_num(optarg);
		else if (c == 'M') mask_tad = 1;
		else if (c == 't') is_tad = 1;
		else if (c == 'g') is_graph = 1;
		else if (c == 'D') is_dedup = 0;
		else if (c == 'v') hk_verbose = atoi(optarg);
	}
	if (argc - optind == 0) {
		print_usage(stderr, &opt);
		return 1;
	}

	m = hk_map_read(argv[optind]);

	if (is_graph) { // for testing only
		g = hk_graph_gen(m, &opt);
		if (g) hk_graph_destroy(g);
	}
	if (is_pairs || is_tad) {
		if (m->pairs == 0 && m->seg)
			m->pairs = hk_map2pairs(m, &m->n_pairs, opt.min_dist, opt.max_seg, opt.min_mapq);
		if (is_dedup)
			m->n_pairs = hk_pair_dedup(m->n_pairs, m->pairs, opt.min_dist);
		if (is_tad || mask_tad) {
			int32_t n_tads;
			struct hk_pair *tads;
			tads = hk_tad_call(m->d, m->n_pairs, m->pairs, opt.max_radius, opt.area_weight, &n_tads);
			if (mask_tad) {
				m->n_pairs = hk_tad_mask(n_tads, tads, m->n_pairs, m->pairs);
				hk_pair_print(stdout, m->d, m->n_pairs, m->pairs);
			} else hk_pair_print(stdout, m->d, n_tads, tads);
			free(tads);
		} else hk_pair_print(stdout, m->d, m->n_pairs, m->pairs);
	} else if (m->seg) {
		hk_map_print(stdout, m);
	} else if (hk_verbose >= 1) {
		fprintf(stderr, "ERROR: the input is not in the seg format\n");
		return 1;
	}
	hk_map_destroy(m);
	return 0;
}
