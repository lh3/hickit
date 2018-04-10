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

static void print_usage(FILE *fp)
{
	fprintf(fp, "Usage: hickit [options] <in.seg>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -p       output pairs\n");
	fprintf(fp, "  -D       retain dup pairs\n");
	fprintf(fp, "  -r NUM   max radius [10m]\n");
}

int main(int argc, char *argv[])
{
	struct hk_opt opt;
	struct hk_map *m = 0;
	struct hk_graph *g = 0;
	int c, is_pairs = 0, is_graph = 0, is_dedup = 1, is_tad = 0;

	hk_opt_init(&opt);
	while ((c = getopt(argc, argv, "pgtDr:v:d:")) >= 0) {
		if (c == 'p') is_pairs = 1;
		else if (c == 't') is_tad = 1;
		else if (c == 'g') is_graph = 1;
		else if (c == 'D') is_dedup = 0;
		else if (c == 'r') opt.max_radius = hk_parse_num(optarg);
		else if (c == 'd') opt.min_dist = hk_parse_num(optarg);
		else if (c == 'v') hk_verbose = atoi(optarg);
	}
	if (argc - optind == 0) {
		print_usage(stderr);
		return 1;
	}

	m = hk_map_read(argv[optind]);

	if (is_graph) { // for testing only
		g = hk_graph_gen(m, &opt);
		if (g) hk_graph_destroy(g);
	}
	if (is_pairs || is_tad) {
		int32_t n_pairs;
		struct hk_pair *pairs;
		pairs = hk_map2pairs(m, &n_pairs, opt.min_dist, opt.max_seg, opt.min_mapq);
		if (is_dedup)
			n_pairs = hk_pair_dedup(n_pairs, pairs, opt.min_dist);
		if (is_tad) {
			hk_tad_call(m->d, n_pairs, pairs, opt.max_radius, 0.0f, NULL);
		} else hk_pair_print(stdout, m->d, n_pairs, pairs);
		free(pairs);
	} else {
		hk_map_print(stdout, m);
	}
	hk_map_destroy(m);
	return 0;
}
