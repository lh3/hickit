#include <stdio.h>
#include <unistd.h>
#include "hickit.h"

int main_view(int argc, char *argv[])
{
	struct hk_map *m;
	int c, is_pairs = 0;

	while ((c = getopt(argc, argv, "p")) >= 0) {
		if (c == 'p') is_pairs = 1;
	}
	if (argc - optind == 0) {
		fprintf(stderr, "Usage: hickit view [options] <in.seg>\n");
		return 1;
	}

	m = hk_map_read(argv[optind]);
	if (is_pairs) hk_map_print_pairs(stdout, m, -1, 3, 20);
	else hk_map_print(stdout, m);
	hk_map_destroy(m);
	return 0;
}

int main(int argc, char *argv[])
{
	return main_view(argc, argv);
}
