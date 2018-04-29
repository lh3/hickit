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

#define hk_ppos1(p) ((int32_t)((p)->pos>>32))
#define hk_ppos2(p) ((int32_t)(p)->pos)
#define hk_intra(p) ((int32_t)((p)->chr>>32) == (int32_t)((p)->chr))

#ifdef __cplusplus
extern "C" {
#endif

int hk_pair_is_sorted(int32_t n_pairs, const struct hk_pair *pairs);
void hk_pair_sort(int32_t n_pairs, struct hk_pair *pairs);

#ifdef __cplusplus
}
#endif

#endif
