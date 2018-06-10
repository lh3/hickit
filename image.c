#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "hkpriv.h"
#include "ksort.h"
KSORT_INIT_GENERIC(float)

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

struct cnt_aux {
	double cnt[4];
	uint32_t tot;
};

void hk_pair_image(const struct hk_sdict *d, int32_t n_pairs, const struct hk_pair *pairs, int w, float phase_thres, int no_grad, const char *fn)
{
	int64_t tot_len, *off;
	int32_t i, j, ww = w * w, pixel_bp, m_tmp, n_tmp;
	uint8_t *buf;
	float *tmp, c_max;
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
	for (i = 0; i < n_pairs; ++i) {
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
		for (j = 0; j < 4; ++j)
			c->cnt[j] += p->_.p4[j];
	}

	// figure out the max depth
	m_tmp = n_tmp = 0, tmp = 0;
	for (i = 0; i < w; ++i) {
		for (j = i; j < w; ++j) {
			int32_t z = i * w + j;
			struct cnt_aux *c = &cnt[z];
			if (c->tot == 0) continue;
			if (phase_thres > 0.0) {
				int32_t k;
				double max = -1.0;
				for (k = 0; k < 4; ++k)
					max = max > c->cnt[k]? max : c->cnt[k];
				if (max / c->tot >= phase_thres) {
					if (n_tmp == m_tmp)
						EXPAND(tmp, m_tmp);
					tmp[n_tmp++] = max;
				}
			} else {
				if (n_tmp == m_tmp)
					EXPAND(tmp, m_tmp);
				tmp[n_tmp++] = c->tot;
			}
		}
	}
	c_max = ks_ksmall_float(n_tmp, tmp, (int)(n_tmp * .95 + .499));
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
				if (phase_thres > 0.0) {
					int x, k, max_k = -1;
					double max = -1.0;
					for (k = 0; k < 4; ++k)
						if (c->cnt[k] > max)
							max_k = k, max = c->cnt[k];
					x = max >= c->tot * phase_thres? (int)(max * t) : (int)(c->tot * t);
					if (no_grad || x > 255) x = 255;
					if (max < c->tot * phase_thres) *p++ = x/2, *p++ = x/2, *p++ = x/2, *q++ = x/2, *q++ = x/2, *q++ = x/2;
					else if (max_k == 0) *p++ = x, *p++ = 0, *p++ = 0, *q++ = x, *q++ = 0, *q++ = 0;
					else if (max_k == 1) *p++ = x, *p++ = 0, *p++ = x, *q++ = 0, *q++ = x, *q++ = x;
					else if (max_k == 2) *p++ = 0, *p++ = x, *p++ = x, *q++ = x, *q++ = 0, *q++ = x;
					else if (max_k == 3) *p++ = 0, *p++ = x, *p++ = 0, *q++ = 0, *q++ = x, *q++ = 0;
				} else {
					int x = (int)(c->tot * t);
					if (no_grad || x > 255) x = 255;
					*p++ = x, *p++ = x, *p++ = x, *q++ = x, *q++ = x, *q++ = x;
				}
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
