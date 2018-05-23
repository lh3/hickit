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

struct hk_sdict *hk_sd_init(void);
int32_t hk_sd_get(const struct hk_sdict *d, const char *s);
int32_t hk_sd_put(struct hk_sdict *d, const char *s, int32_t len);
struct hk_sdict *hk_sd_split_phase(const struct hk_sdict *d, int32_t *ploidy_XY);
void hk_sd_destroy(struct hk_sdict *d);

int hk_pair_is_sorted(int32_t n_pairs, const struct hk_pair *pairs);
void hk_pair_sort(int32_t n_pairs, struct hk_pair *pairs);

void hk_bmap_set_offcnt(struct hk_bmap *m);
int hk_bmap_pos2bid(const struct hk_bmap *m, int32_t chr, int32_t pos);

int32_t ks_ksmall_int32_t(size_t n, int32_t arr[], size_t kk);
float ks_ksmall_float(size_t n, float arr[], size_t kk);

static inline uint64_t hash64(uint64_t key)
{
	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return key;
}

#ifdef __cplusplus
}
#endif

#endif
