/* The MIT License

   Copyright (c) 2008-2009, by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef _AC_KLIST_H
#define _AC_KLIST_H

#include <stdlib.h>

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */

#define KMEMPOOL_INIT2(SCOPE, name, kmptype_t, kmpfree_f)				\
	typedef struct {													\
		size_t cnt, n, max;												\
		kmptype_t **buf;												\
	} kmp_##name##_t;													\
	SCOPE kmp_##name##_t *kmp_init_##name(void) {						\
		return (kmp_##name##_t*)calloc(1, sizeof(kmp_##name##_t));		\
	}																	\
	SCOPE void kmp_destroy_##name(kmp_##name##_t *mp) {					\
		size_t k;														\
		for (k = 0; k < mp->n; ++k) {									\
			kmpfree_f(mp->buf[k]); free(mp->buf[k]);					\
		}																\
		free(mp->buf); free(mp);										\
	}																	\
	SCOPE kmptype_t *kmp_alloc_##name(kmp_##name##_t *mp) {				\
		++mp->cnt;														\
		if (mp->n == 0) return (kmptype_t*)calloc(1, sizeof(kmptype_t)); \
		return mp->buf[--mp->n];										\
	}																	\
	SCOPE void kmp_free_##name(kmp_##name##_t *mp, kmptype_t *p) {		\
		--mp->cnt;														\
		if (mp->n == mp->max) {											\
			mp->max = mp->max? mp->max<<1 : 16;							\
			mp->buf = (kmptype_t**)realloc(mp->buf, sizeof(kmptype_t *) * mp->max); \
		}																\
		mp->buf[mp->n++] = p;											\
	}

#define KMEMPOOL_INIT(name, kmptype_t, kmpfree_f)						\
	KMEMPOOL_INIT2(static inline klib_unused, name, kmptype_t, kmpfree_f)

#define kmempool_t(name) kmp_##name##_t
#define kmp_init(name) kmp_init_##name()
#define kmp_destroy(name, mp) kmp_destroy_##name(mp)
#define kmp_alloc(name, mp) kmp_alloc_##name(mp)
#define kmp_free(name, mp, p) kmp_free_##name(mp, p)

#endif
