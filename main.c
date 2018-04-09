#include <stdio.h>
#include "hickit.h"

int main(int argc, char *argv[])
{
	struct hk_map *m;
	if (argc == 1) {
		return 1;
	}
	m = hk_map_read(argv[1]);
	hk_map_print(stdout, m);
	hk_map_destroy(m);
	return 0;
}
