#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
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

int main(int argc, char *argv[])
{
	struct hk_gopt opt;
	struct hk_map *m = 0;
	struct hk_graph *g = 0;
	int c, is_pairs = 0, is_graph = 0;

	hk_gopt_init(&opt);
	while ((c = getopt(argc, argv, "pgr:")) >= 0) {
		if (c == 'p') is_pairs = 1;
		else if (c == 'g') is_graph = 1;
		else if (c == 'r') opt.max_radius = hk_parse_num(optarg);
	}
	if (argc - optind == 0) {
		fprintf(stderr, "Usage: hickit [options] <in.seg>\n");
		return 1;
	}

	m = hk_map_read(argv[optind]);
	if (is_graph) {
		g = hk_graph_gen(m, &opt);
	} else if (is_pairs) {
		hk_map_print_pairs(stdout, m, -1, 3, 20);
	} else {
		hk_map_print(stdout, m);
	}
	if (g) hk_graph_destroy(g);
	hk_map_destroy(m);
	return 0;
}
