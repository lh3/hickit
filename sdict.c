#include <string.h>
#include "hkpriv.h"
#include "khash.h"
KHASH_MAP_INIT_STR(str, int32_t)
typedef khash_t(str) sdict_t;

struct hk_sdict *hk_sd_init(void)
{
	struct hk_sdict *d;
	d = CALLOC(struct hk_sdict, 1);
	d->h = kh_init(str);
	return d;
}

void hk_sd_destroy(struct hk_sdict *d)
{
	int32_t i;
	for (i = 0; i < d->n; ++i)
		free(d->name[i]);
	free(d->name);
	free(d->len);
	kh_destroy(str, (sdict_t*)d->h);
}

int32_t hk_sd_put(struct hk_sdict *d, const char *s, int32_t len)
{
	sdict_t *h = (sdict_t*)d->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, s, &absent);
	if (absent) {
		if (d->n == d->m) {
			EXPAND(d->name, d->m);
			REALLOC(d->len, d->m);
		}
		d->len[d->n] = len;
		kh_key(h, k) = d->name[d->n] = strdup(s);
		kh_val(h, k) = d->n++;
	}
	return kh_val(h, k);
}

struct hk_sdict *hk_sd_sep_phase(const struct hk_sdict *d, int32_t *ploidy_XY)
{
	struct hk_sdict *dp;
	int i;
	dp = hk_sd_init();
	for (i = 0; i < d->n; ++i) {
		int ploidy = ploidy_XY? ploidy_XY[i] >> 8 : 1;
		if (ploidy <= 1) {
			hk_sd_put(dp, d->name[i], d->len[i]);
		} else {
			int l, j;
			char *s;
			l = strlen(d->name[i]);
			s = CALLOC(char, l + 2);
			strcpy(s, d->name[i]);
			s[l+1] = 0;
			for (j = 0; j < ploidy; ++j) {
				s[l] = 'a' + j;
				hk_sd_put(dp, s, d->len[i]);
			}
			free(s);
		}
	}
	return dp;
}

struct hk_sdict *hk_sd_dup(const struct hk_sdict *d)
{
	return hk_sd_sep_phase(d, 0);
}

int32_t *hk_sd_ploidy_XY(const struct hk_sdict *d, int32_t *sex_flag)
{
	int32_t i, *ploidy, n_X = 0, n_Y = 0;
	if (sex_flag) *sex_flag = 0;
	for (i = 0; i < d->n; ++i) {
		int32_t is_X, is_Y;
		is_X = (strchr(d->name[i], 'X') != 0);
		is_Y = (strchr(d->name[i], 'Y') != 0);
		if (is_X && is_Y) continue;
		if (is_X) ++n_X;
		if (is_Y) ++n_Y;
	}
	if (n_X == 0 && n_Y == 0)
		fprintf(stderr, "[W::%s] no sex chromosomes (identified by 'X' or 'Y' in chr names)\n", __func__);
	if (n_X > 1 || n_Y > 1) {
		fprintf(stderr, "[E::%s] multiple chr contain 'X' or 'Y' in names\n", __func__);
		exit(1);
	}
	if (n_X > 0 && sex_flag) *sex_flag |= 1;
	if (n_Y > 0 && sex_flag) *sex_flag |= 2;
	ploidy = CALLOC(int32_t, d->n);
	for (i = 0; i < d->n; ++i) {
		int32_t is_X, is_Y;
		is_X = (strchr(d->name[i], 'X') != 0);
		is_Y = (strchr(d->name[i], 'Y') != 0);
		if (is_X && is_Y) continue;
		if (n_Y == 0) ploidy[i] = 2 << 8;
		else ploidy[i] = (is_X || is_Y? 1 : 2) << 8;
		if (is_X) ploidy[i] |= 1;
		if (is_Y) ploidy[i] |= 2;
	}
	return ploidy;
}

