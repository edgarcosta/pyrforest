#ifndef _ZZ_MEM_INCLUDE_
#define _ZZ_MEM_INCLUDE_

#include <stdlib.h>

#define ZZ_MEM_TRACKING     1
#define ZZ_MEM_THREADING    1
#define ZZ_MEM_LOGGING      0

extern long zz_mem_used;
extern long zz_mem_peak;
extern long zz_malloc_count;
extern long zz_free_count;

#if ZZ_MEM_TRACKING
void *zz_malloc (size_t n);
void zz_free (void *p, size_t n);
#else
static inline void *zz_malloc (size_t n) { return malloc (n); }
static inline void zz_free (void *p, size_t n) { free (p); }
#endif

#endif
