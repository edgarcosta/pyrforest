#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <memory.h>
#include "hwmem.h"
#include "zzmem.h"

#if HW_MEM_TRACKING
long hw_mem_used;
long hw_mem_peak;
long hw_malloc_count;
long hw_realloc_count;
long hw_free_count;
long zz_overhead;

void hw_mem_init (size_t n)
{
}

void hw_mem_clear (void)
{
}

void *hw_malloc (size_t n)
{
    void *p;
#if HW_MEM_THREADING
    #pragma omp critical
#endif
{
    hw_mem_used += n;
    if ( hw_mem_used > hw_mem_peak ) hw_mem_peak = hw_mem_used;
    hw_malloc_count++;
}
    p = malloc (n);
#if HW_MEM_LOG
    printf ("hw_malloc %ld (%lx)\n", n, (unsigned long)p);
#endif
    return p;
}
void *hw_calloc (size_t n)
{
    void *p;
#if HW_MEM_THREADING
    #pragma omp critical
#endif
{
    hw_mem_used += n;
    if ( hw_mem_used > hw_mem_peak ) hw_mem_peak = hw_mem_used;
    hw_malloc_count++;
}
    p = calloc (n,1);
#if HW_MEM_LOG
    printf ("hw_calloc %ld (%lx)\n", n, (unsigned long)p);
#endif
    return p;
}

void *hw_realloc (void *p, size_t n1, size_t n2)
{
    void *newp;
#if HW_MEM_THREADING
    #pragma omp critical
#endif
{
    hw_mem_used += n2;
    hw_mem_used -= n1;
    if ( hw_mem_used > hw_mem_peak ) hw_mem_peak = hw_mem_used;
    hw_realloc_count++;
}
    newp = realloc (p, n2);
#if HW_MEM_LOG
printf ("realloc %ld (%lx) to %ld (%lx)\n", n1, (unsigned long)p, n2, (unsigned long)newp);
#endif
    return newp;
}

void hw_free (void *p, size_t n)
{
#if HW_MEM_LOG
printf ("hw_free %ld (%lx)\n", n, (unsigned long)p);
#endif
#if HW_MEM_THREADING
    #pragma omp critical
#endif
{
    hw_mem_used -= n;
    hw_free_count++;
}
    free (p);
}

double hw_peak (void)
{
    return (hw_mem_peak+(zz_mem_peak-zz_overhead))/(1<<20);
}

void hw_mem_report (int verbosity)
{
    if ( verbosity < 0 ) return;
#if HW_MEM_THREADING
    #pragma omp critical
#endif
{
    if ( !verbosity && !hw_mem_used && hw_malloc_count == hw_free_count && !zz_mem_used) {
        printf ("%.3f MB memory peak usage (of which %.3f MB is in zz, excluding %.3f MB fixed overhead)\n", (double)(hw_mem_peak+(zz_mem_peak-zz_overhead))/(1<<20), (double)(zz_mem_peak-zz_overhead)/(1<<20), (double)zz_overhead/(1<<20));
    } else {
        printf ("hw_mem_report\n");
        printf ("    %ld allocs, %ld reallocs, %ld frees\n",
            hw_malloc_count, hw_realloc_count, hw_free_count);
        printf ("    %ld bytes outstanding in %ld blocks\n", hw_mem_used,
            hw_malloc_count - hw_free_count);
        printf ("    %ld bytes peak usage (excluding heap)\n", hw_mem_peak);
        printf ("zz_mem_report\n");
        printf ("    %ld allocs, %ld frees\n",
            zz_malloc_count, zz_free_count);
        printf ("    %ld bytes outstanding in %ld blocks\n", zz_mem_used,
            zz_malloc_count - zz_free_count);
        printf ("    %ld bytes peak usage beyond %ld fixed overhead\n", zz_mem_peak-zz_overhead, zz_overhead);
        printf ("%.3f MB memory peak usage (of which %.3f MB is in zz, excluding %.3f MB fixed overhead)\n", (double)(hw_mem_peak+(zz_mem_peak-zz_overhead))/(1<<20), (double)(zz_mem_peak-zz_overhead)/(1<<20), (double)zz_overhead/(1<<20));
    }
}
}
#endif
