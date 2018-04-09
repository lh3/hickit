#ifndef HICKIT_H
#define HICKIT_H

#include <stdint.h>

#define HK_SUB_DELIM '!'

#ifdef __cplusplus
extern "C" {
#endif

struct hk_sdict {
	int32_t n, m;
	char **name;
	void *h;
};

struct hk_seg {
	int32_t frag_id, chr, st, en;
	int8_t strand, phase;
	int16_t mapq;
};

struct hk_map {
	struct hk_sdict *d;
	int32_t n_frag, m_seg, n_seg;
	struct hk_seg *seg;
};

struct hk_pair {
	uint64_t pos[2];
	int8_t phase[2];
	int8_t rel_strand;
};

struct hk_map *hk_map_read(const char *fn);
void hk_map_destroy(struct hk_map *m);
void hk_map_print(FILE *fp, const struct hk_map *m);

#ifdef __cplusplus
}
#endif

#endif
