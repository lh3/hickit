#ifndef HK_PRIV_H
#define HK_PRIV_H

#include <stdlib.h>
#include "hickit.h"

#define MALLOC(type, len) ((type*)malloc((len) * sizeof(type)))
#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#define EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((a), (m)); \
	} while (0)

#ifdef __cplusplus
extern "C" {
#endif

struct hk_pair *hk_map2pairs(const struct hk_map *m, int32_t *_n_pairs, int min_dist, int max_seg, int min_mapq);
int32_t hk_pair_dedup(int n_pairs, struct hk_pair *pairs);

#ifdef __cplusplus
}
#endif

#endif
