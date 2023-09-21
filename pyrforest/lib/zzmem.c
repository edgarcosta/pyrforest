#include "zzmem.h"

long zz_mem_used;
long zz_mem_peak;
long zz_malloc_count;
long zz_free_count;

#if ZZ_MEM_TRACKING

void *zz_malloc (size_t n)
{
	void *p;
#if ZZ_MEM_THREADING
	#pragma omp critical
#endif
	{
		zz_mem_used += n;
		if ( zz_mem_used > zz_mem_peak ) zz_mem_peak = zz_mem_used;
		zz_malloc_count++;
	}
	p = malloc (n);
#if ZZ_MEM_LOGGING
	printf ("zz_malloc %ld (%lx)\n", n, (unsigned long)p);
#endif
	return p;
}


void zz_free (void *p, size_t n)
{
#if ZZ_MEM_LOGGING
printf ("zz_free %ld (%lx)\n", n, (unsigned long)p);
#endif
#if ZZ_MEM_THREADING
	#pragma omp critical
#endif
	{
	zz_mem_used -= n;
	zz_free_count++;
	}
	free (p);
}
#endif
