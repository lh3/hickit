#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "hkpriv.h"
#include "ksort.h"
KSORT_INIT_GENERIC(int)

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void hk_pair_image(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int w, const char *fn)
{
	int64_t tot_len, *off;
	int32_t i, *cnt, j, n1, *tmp, ww = w * w, c_max;
	uint8_t *p, *buf;
	double s;
	off = CALLOC(int64_t, d->n);
	for (tot_len = 0, i = 0; i < d->n; ++i) {
		off[i] = tot_len;
		tot_len += d->len[i];
	}
	s = (double)w / tot_len;
	cnt = CALLOC(int32_t, w * w);
	for (i = 0, n1 = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		int64_t x[2], z;
		int32_t y[2];
		x[0] = off[p->chr>>32] + hk_ppos1(p);
		x[1] = off[(int32_t)p->chr] + hk_ppos2(p);
		y[0] = (int32_t)(x[0] * s);
		y[1] = (int32_t)(x[1] * s);
		assert(y[0] < w && y[1] < w);
		z = y[0] * w + y[1];
		++cnt[z];
		if (cnt[z] == 1) ++n1;
	}
	tmp = CALLOC(int32_t, n1);
	for (i = j = 0; i < ww; ++i)
		if (cnt[i] > 0) tmp[j++] = cnt[i];
	assert(j == n1);
	c_max = ks_ksmall_int(n1, tmp, (int)(n1 * .95 + .499));
	free(tmp);
	s = (double)255 / c_max;
	buf = CALLOC(uint8_t, ww * 3);
	for (i = 0, p = buf; i < w; ++i) {
		for (j = 0; j < w; ++j) {
			int32_t z = i * w + j;
			/*
			if (cnt[z] == 0) {
				*p++ = 255, *p++ = 255, *p++ = 255;
			} else {
				int x = (int)(s * cnt[z]);
				if (x > 255) x = 255;
				*p++ = 255, *p++ = 255 - x, *p++ = 255 - x;
			}
			*/
			if (cnt[z] == 0) {
				*p++ = 0, *p++ = 0, *p++ = 0;
			} else {
				int x = (int)(s * cnt[z]);
				if (x > 255) x = 255;
				*p++ = x, *p++ = 0, *p++ = 0;
			}
		}
	}
	stbi_write_png(fn, w, w, 3, buf, w * 3);
	free(buf);
	free(cnt);
	free(off);
}
