#ifndef _HWMEM_INCLUDE_
#define _HWMEM_INCLUDE_

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HW_MEM_LOG          0       // will log every allocation and free, use only for debugging
#define HW_MEM_TRACKING     1       // keeps track of # and total size of all memory allocations
#define HW_MEM_THREADING    1       // make memory tracking thread safe (use only for debugging)

extern long hw_mem_used, hw_mem_peak;

// wrappers to track memory usage
#if HW_MEM_TRACKING
void hw_mem_init (size_t n);
void hw_mem_clear (void);
void *hw_malloc (size_t n);
void *hw_calloc (size_t n);
void *hw_realloc (void *p, size_t n1, size_t n2);
void hw_free (void *p, size_t n);
double hw_peak (void);
void hw_mem_report (int verbosity);

// be sure to call this first thing if you want to track gmp memory usage
static inline void hw_mem_track_gmp (void) { mp_set_memory_functions (hw_malloc, hw_realloc, hw_free); }
#else
static inline void hw_mem_init (size_t n) {}
static inline void hw_mem_clear (void) {}
static inline void *hw_malloc (size_t n) { return malloc (n); }
static inline void *hw_calloc (size_t n) { return calloc (n,1); }
static inline void *hw_realloc (void *p, size_t n1, size_t n2) { return realloc (p, n2); }
static inline void hw_free (void *p, size_t n) { free (p); }
static inline double hw_peak (void) { return 0.0; }
static inline void hw_mem_report (int verbosity) { if ( verbosity ) puts ("HW_MEM_TRACKING = 0"); }
static inline void hw_mem_track_gmp (void) {return; }
#endif

#ifdef __cplusplus
}
#endif

#endif
