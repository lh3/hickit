#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "hkpriv.h"
#include "ksort.h"
KSORT_INIT_GENERIC(int)

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

struct cnt_aux {
	int32_t cnt[4];
	uint32_t tot;
};

void hk_pair_image(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int w, float phase_thres, const char *fn)
{
	int64_t tot_len, *off;
	int32_t i, j, n1, *tmp, ww = w * w, c_max, pixel_bp;
	uint8_t *buf;
	double s, t;
	struct cnt_aux *cnt;

	for (tot_len = 0, i = 0; i < d->n; ++i)
		tot_len += d->len[i];
	s = (double)(w - d->n + 1) / tot_len;
	pixel_bp = tot_len / (w - d->n + 1) + 1;

	off = CALLOC(int64_t, d->n);
	for (tot_len = 0, i = 0; i < d->n; ++i) {
		off[i] = tot_len;
		tot_len += d->len[i] + pixel_bp;
	}

	cnt = CALLOC(struct cnt_aux, w * w);
	for (i = 0, n1 = 0; i < n_pairs; ++i) {
		const struct hk_pair *p = &pairs[i];
		int64_t x[2], z;
		int32_t y[2];
		struct cnt_aux *c;
		x[0] = off[p->chr>>32] + hk_ppos1(p);
		x[1] = off[(int32_t)p->chr] + hk_ppos2(p);
		y[0] = (int32_t)(x[0] * s);
		y[1] = (int32_t)(x[1] * s);
		assert(y[0] < w && y[1] < w);
		z = y[0] * w + y[1];
		c = &cnt[z];
		++c->tot;
		if (phase_thres > 0.0f) {
			int q[2];
			q[0] = p->_.phase_prob[0] <= phase_thres? 0 : p->_.phase_prob[0] >= 1.0f - phase_thres? 1 : -1;
			q[1] = p->_.phase_prob[1] <= phase_thres? 0 : p->_.phase_prob[1] >= 1.0f - phase_thres? 1 : -1;
			if (q[0] >= 0 && q[1] >= 0)
				++c->cnt[q[0]<<1|q[1]];
		}
		if (c->tot == 1) ++n1;
	}

	tmp = CALLOC(int32_t, n1);
	for (i = j = 0; i < ww; ++i)
		if (cnt[i].tot) tmp[j++] = cnt[i].tot;
	assert(j == n1);
	c_max = ks_ksmall_int(n1, tmp, (int)(n1 * .99 + .499));
	free(tmp);

	t = 255.0 / c_max;
	buf = CALLOC(uint8_t, ww * 3);
	for (i = 0; i < w; ++i) {
		for (j = i; j < w; ++j) {
			int32_t z = i * w + j;
			uint8_t *q = &buf[(j * w + i) * 3];
			uint8_t *p = &buf[(i * w + j) * 3];
			struct cnt_aux *c = &cnt[z];
			if (c->tot == 0) {
				*p++ = 0, *p++ = 0, *p++ = 0, *q++ = 0, *q++ = 0, *q++ = 0;
			} else {
				int x, k, max_k = -1, max = -1;
				for (k = 0; k < 4; ++k)
					if (c->cnt[k] > max)
						max_k = k, max = c->cnt[k];
				if (max == 0 || max < c->tot/2) x = (int)(c->tot * t);
				else x = (int)(max * t);
				if (x > 255) x = 255;
				if (max == 0 || max < c->tot/2) *p++ = x/2, *p++ = x/2, *p++ = x/2, *q++ = x/2, *q++ = x/2, *q++ = x/2;
				else if (max_k == 0) *p++ = x, *p++ = 0, *p++ = 0, *q++ = x, *q++ = 0, *q++ = 0;
				else if (max_k == 1) *p++ = x, *p++ = 0, *p++ = x, *q++ = 0, *q++ = x, *q++ = x;
				else if (max_k == 2) *p++ = 0, *p++ = x, *p++ = x, *q++ = x, *q++ = 0, *q++ = x;
				else if (max_k == 3) *p++ = 0, *p++ = x, *p++ = 0, *q++ = 0, *q++ = x, *q++ = 0;
			}
		}
	}

	// write lines
	for (i = 0; i < d->n - 1; ++i) {
		int x = (int)(off[i+1] * s + .499);
		uint8_t *p;
		for (j = 0, p = &buf[x * w * 3]; j < w; ++j, p += 3)
			p[0] = 32, p[1] = 32, p[2] = 32;
		for (j = 0, p = &buf[x * 3]; j < w; ++j, p += w * 3)
			p[0] = 32, p[1] = 32, p[2] = 32;
	}

	stbi_write_png(fn, w, w, 3, buf, w * 3);

	free(buf);
	free(cnt);
	free(off);
}
