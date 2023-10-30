#ifndef _HWMPZ_INCLUDE_
#define _HWMPZ_INCLUDE_

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>
#include "hwmem.h"
#include "hwmpz_tune.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int hw_disable_fft;

static inline int _ui_len (uint64_t x) { return x?64-__builtin_clzl(x):0; }
static inline int _si_len (int64_t x) { return 64-__builtin_clzl(x<0?-x:x); }

void hw_mpz_setup ();
void hw_mpz_clear ();

// handy mpz inlines
static inline int mpz_is_zero (mpz_t x) { return mpz_sgn(x) ? 0 : 1; }
static inline void mpz_set_zero (mpz_t o) { mpz_set_ui (o,0); }
static inline int mpz_is_one (mpz_t x) { return mpz_cmp_ui(x,1) == 0 ? 1 : 0; }
static inline void mpz_set_one (mpz_t o) { mpz_set_ui (o,1); }
static inline void mpz_set_negone (mpz_t o) { mpz_set_si (o,-1); }
static inline void mpz_init_one (mpz_t o) { mpz_init_set_ui (o,1); }
static inline int mpz_equal (mpz_t a, mpz_t b) { return mpz_cmp(a,b) == 0; }
static inline long mpz_bits (mpz_t a) { return mpz_sizeinbase(a,2); }
static inline int mpz_remove_ui (mpz_t z, uint64_t f) { int n; for ( n = 0 ; mpz_divisible_ui_p(z,f) ; n++ ) { mpz_divexact_ui(z,z,f); } return n; }

// vector versions of basic mpz_t functions
static inline mpz_t *mpz_vec_alloc (long n) { return (mpz_t*)hw_malloc (n*sizeof(mpz_t)); }
static inline mpz_t *mpz_vec_init (mpz_t *A, long n) { long i; for ( i = 0 ; i < n ; i++ ) { mpz_init (A[i]); } return A; }
static inline mpz_t *mpz_vec_init2 (mpz_t *A, long n, int bits) { long i; for ( i = 0 ; i < n ; i++ ) { mpz_init2 (A[i], bits); } return A; }
static inline mpz_t *mpz_vec_alloc_and_init (long n) { mpz_t *A = (mpz_t*)hw_malloc (n*sizeof(mpz_t)); mpz_vec_init (A,n); return A; }
static inline mpz_t *mpz_vec_alloc_and_init2 (long n, long bits) { mpz_t *A = (mpz_t*)hw_malloc (n*sizeof(mpz_t)); mpz_vec_init2 (A,n,bits); return A; }
static inline mpz_t *mpz_vec_set (mpz_t *A, mpz_t *B, long n) { if ( A == B ) return A; for ( long i = 0 ; i < n ; i++ ) { mpz_set (A[i], B[i]); } return A; }
static inline mpz_t *mpz_vec_set_si (mpz_t *A, long *a, long n) { for ( long i = 0 ; i < n ; i++ ) { mpz_set_si (A[i], a[i]); } return A; }
static inline mpz_t *mpz_vec_set_ui (mpz_t *A, unsigned long *a, long n) { for ( long i = 0 ; i < n ; i++ ) { mpz_set_ui (A[i], a[i]); } return A; }
static inline mpz_t *mpz_vec_set_zero (mpz_t *A, long n) { long i; for ( i = 0 ; i < n ; i++ ) { mpz_set_zero (A[i]); } return A; }
static inline int mpz_vec_equal (mpz_t *A, mpz_t *B, long n) { if ( A == B ) return 1; for ( long i = 0 ; i < n ; i++ ) if ( ! mpz_equal (A[i], B[i]) ) return 0;  return 1; }
static inline mpz_t *mpz_vec_addto (mpz_t *A, mpz_t *B, long n) { for ( long i = 0 ; i < n ; i++ ) { mpz_add (A[i], A[i], B[i]); } return A; }
static inline void mpz_vec_dot_product (mpz_t z, mpz_t *A, mpz_t *B, long n) { mpz_set_zero (z); for ( long i = 0 ; i < n ; i++ ) mpz_addmul (z, A[i], B[i]); }
static inline mpz_t *mpz_vec_scalar_product (mpz_t *A, mpz_t b, mpz_t *B, long n) { for ( long i = 0 ; i < n ; i++ ) { mpz_mul (A[i], b, B[i]); } return A; }
    
// max size of vector entry
static inline long mpz_vec_max_size (mpz_t *A, long n) { long s, t;  s = 0; for ( long i = 0 ; i < n ; i++ ) { t = mpz_size(A[i]); if ( t > s ) s = t; } return s; }
static inline long mpz_vec_total_size (mpz_t *A, long n) { long s;  s = 0; for ( long i = 0 ; i < n ; i++ ) s += mpz_size(A[i]);   return s;}

static inline mpz_t *mpz_vec_mod_naive (mpz_t *A, mpz_t *B, long n, mpz_t m) { for (  long i = 0 ; i < n ; i++ ) { mpz_fdiv_r (A[i], B[i], m); } return A; }
mpz_t *mpz_vec_mod_fft (mpz_t *A, mpz_t *B, long n, mpz_t m, int *reps);
static inline mpz_t *mpz_vec_mod_init_naive (mpz_t *A, mpz_t *B, long n, mpz_t m, mpz_t w)
    { for (  long i = 0 ; i < n ; i++ ) { mpz_fdiv_r (w, B[i], m); mpz_init_set (A[i], w); }  return A; }
mpz_t *mpz_vec_mod_init_fft (mpz_t *A, mpz_t *B, long n, mpz_t m, mpz_t w);
static inline mpz_t *mpz_vec_div_mod_naive (mpz_t *Q, mpz_t *R, mpz_t *A, long n, mpz_t m) { for (  long i = 0 ; i < n ; i++ ) { mpz_fdiv_qr (Q[i], R[i], A[i], m); } return Q; }

static inline void mpz_vec_clear (mpz_t *A, long n) {  int i; for ( i = 0 ; i < n ; i++ ) mpz_clear (A[i]); }
static inline void mpz_vec_free (mpz_t *A, long n) { hw_free (A,n*sizeof(mpz_t) ); }
static inline void mpz_vec_clear_and_free (mpz_t *A, long n) { mpz_vec_clear (A, n);  mpz_vec_free (A, n); }

// computes x = c*\prod_{i=0}^d fpow[i][fdeg[i]]
static inline void mpz_monomial_eval_pow (mpz_t x, long c, unsigned char *fdeg, mpz_t **fpow, int d)
{
    int i;
    
    if ( c < 0 ) { mpz_set_ui (x,-c); mpz_neg (x,x); } else { mpz_set_ui(x,c); }
    for ( i = 0 ; i <= d ; i++ ) if ( fdeg[i] ) mpz_mul (x,x,fpow[i][fdeg[i]]);
}

// computes y = g(x), given an array of powers of x (same as dot product, except xpow[0] is always 1)
static inline void mpz_poly_eval_pow (mpz_t y, mpz_t g[], int d, mpz_t xpow[])
    { mpz_set (y, g[0]);  for (  int i = 1 ; i <= d ; i++ ) mpz_addmul (y, g[i], xpow[i]); }

static inline mpz_t *mpz_matrix_alloc (int d) { return mpz_vec_alloc (d*d); }
static inline void mpz_matrix_free (mpz_t *A, int d) { mpz_vec_free (A, d*d); }
static inline void mpz_matrix_init (mpz_t *M, int d) { mpz_vec_init (M, d*d); }
static inline mpz_t *mpz_matrix_alloc_and_init (int d) { return mpz_vec_alloc_and_init (d*d); }
static inline void mpz_matrix_init2 (mpz_t *M, int d, long b) { mpz_vec_init2 (M, d*d, b); }
static inline mpz_t *mpz_matrix_alloc_and_init2 (mpz_t *M, int d, long b) { return mpz_vec_alloc_and_init2 (d*d, b); }
static inline void mpz_matrix_set_zero (mpz_t *M, int d) { mpz_vec_set_zero (M, d*d); }
static inline void mpz_matrix_set_one (mpz_t *M, int d) { mpz_vec_set_zero (M, d*d); for ( int i = 0 ; i < d ; i++ ) mpz_set_one (M[i*d+i]); }
static inline void mpz_matrix_clear (mpz_t *M, int d) { mpz_vec_clear (M, d*d); }
static inline void mpz_matrix_clear_and_free (mpz_t *M, int d) { mpz_vec_clear_and_free (M, d*d); }
static inline int mpz_matrix_equal (mpz_t *A, mpz_t *B, int d) { return mpz_vec_equal (A, B, d*d); }
static inline void mpz_matrix_set (mpz_t *A, mpz_t *B, int d) { mpz_vec_set (A, B, d*d); }
static inline void mpz_matrix_div (mpz_t *A, int d, mpz_t mdiv) { for (  int i = 0 ; i < d*d ; i++ ) mpz_divexact (A[i], A[i], mdiv); } 

// handy print functions for debugging
static inline void _mpz_row_print (mpz_t *A, int d, int w)
{
    int i;
    
    printf ("[ "); 
    for ( i = 0 ; i < d ; i++ ) gmp_printf ("%*Zd ", w+1, A[i]);
    printf ("]");
}

static inline void mpz_row_print (mpz_t *A, int d)
{
    int i, x, w;
    
    for ( w = 1, i = 0 ; i < d ; i++ ) { x = mpz_sizeinbase(A[i],10);  if ( x > w ) w = x; }
    _mpz_row_print (A, d, w);
}

static inline void _i_row_print (long *A, int d, int w)
{
     int i;
    
    printf ("[ "); 
    for ( i = 0 ; i < d ; i++ ) printf ("%*ld ", w, A[i]);
    printf ("]");
}

static inline void i_row_print (long *A, int d)
{
    char buf[32];
    int i, x, w;
    
    for ( w = 1, i = 0 ; i < d ; i++ ) { x = sprintf(buf,"%ld",A[i]);  if ( x > w ) w = x; }
    _i_row_print (A, d, w);
}

static inline void mpz_row_print_sizes (mpz_t *A, int d)
{
    int i;
    
    printf ("[ "); 
    for ( i = 0 ; i < d ; i++ ) gmp_printf ("%6ld ", mpz_size(A[i]));
    printf (" ]");
}

static inline void mpz_matrix_print (mpz_t *M, int d)
{
    int i,j,w,x;

    for ( w=1, i = 0 ; i < d ; i++ ) for ( j = 0 ; j < d ; j++ ) { x = mpz_sizeinbase(M[i*d+j],10);  if ( x > w ) w = x; }
    for ( i = 0 ; i < d ; i++ ) { _mpz_row_print (M+i*d, d, w); printf ("\n"); }
}

static inline void i_matrix_print (long *M, int d)
{
    char buf[32];
    int i,j,w,x;

    for ( w=1, i = 0 ; i < d ; i++ ) for ( j = 0 ; j < d ; j++ ) { x = sprintf(buf,"+%ld",M[i*d+j]);  if ( x > w ) w = x; }
    for ( i = 0 ; i < d ; i++ ) { _i_row_print (M+i*d, d, w); printf ("\n"); }
}

static inline void mpz_matrix_print_sizes (mpz_t *M, int d)
{
    int i;

    for ( i = 0 ; i < d ; i++ ) { mpz_row_print_sizes (M+i*d, d); puts(""); }
}

// computes C=AB mod m, where A and C are r-by-d matrices and B, is a d-by-d matrices, using naive rd^2 alg, w is a work variable (helps to avoid reallocs during computation)
// ALIASING NOT ALLOWED
static inline void mpz_rmatrix_mult_naive (mpz_t *C, mpz_t *A, int r, mpz_t *B, int d, mpz_t w)
{
    for ( int i = 0 ; i < r ; i++ ) for ( int j = 0 ; j < d ; j++ ) {
        mpz_mul (w,A[i*d],B[j]);
        for (  int k = 1 ; k < d ; k++ ) mpz_addmul (w,A[i*d+k],B[k*d+j]);
        mpz_set (C[i*d+j],w);
    }
}
static inline void mpz_matrix_mult_naive (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t w) { mpz_rmatrix_mult_naive (C, A, d, B, d, w); }

// computes C=AB, where A, B, and C are square d-by-d matrices, using naive d^3 alg, w is a work variable (helps to avoid reallocs during computation)
// initializes the entries of C as it computes them
// NO ALIASING
mpz_t *mpz_rmatrix_mult_fft (mpz_t *C, mpz_t *A, int r, mpz_t *B, int d, mpz_t w, int *reps1, int *reps2);
static inline mpz_t *mpz_row_matrix_mult_fft (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t w) { return mpz_rmatrix_mult_fft (C, A, 1, B, d, w, NULL, NULL); }
static inline mpz_t *mpz_matrix_mult_fft (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t w) { return mpz_rmatrix_mult_fft (C, A, d, B, d, w, NULL, NULL); }

// computes C=AB mod m, where A and C are r-by-d matrices and B, is a d-by-d matrices, using naive rd^2 alg, w is a work variable (helps to avoid reallocs during computation)
// ALIASING NOT ALLOWED
static inline void mpz_rmatrix_mult_mod_naive (mpz_t *C, mpz_t *A, int r, mpz_t *B, int d, mpz_t m, mpz_t w)
{
    int i, j, k;
    
    for ( i = 0 ; i < r ; i++ ) for ( j = 0 ; j < d ; j++ ) {
        mpz_mul (w,A[i*d],B[j]);
        for ( k = 1 ; k < d ; k++ ) mpz_addmul (w,A[i*d+k],B[k*d+j]);
        mpz_fdiv_r (C[i*d+j],w,m);
    }
}
static inline void mpz_matrix_mult_mod_naive (mpz_t *C, mpz_t *A, int r, mpz_t *B, int d, mpz_t m, mpz_t w) { mpz_rmatrix_mult_mod_naive (C, A, r, B, d, m, w); }

// computes C=AB, where B is a square d-by-d matrix, while A and C are d-entry row vectors.
// NO ALIASING
static inline void mpz_row_matrix_mult_naive (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t w)
{
    int i, j;
    
    for ( i = 0 ; i < d ; i++ ) {
        mpz_mul (w,A[0],B[i]);
        for ( j = 1 ; j < d ; j++ ) mpz_addmul (w,A[j],B[j*d+i]);
        mpz_set (C[i],w);
    }
}

// return an upper boiund on the number of bits in x (for efficiency just count limbs)
static inline long mpz_height (mpz_t x) { return GMP_NUMB_BITS*mpz_size(x); }
static inline long mpz_vec_height (mpz_t *A, long n) { return GMP_NUMB_BITS*mpz_vec_max_size(A,n); }

// returns an upper bound on the number of bits in the largest entry of A
static inline long mpz_row_height (mpz_t *A, int d) { return mpz_vec_height (A, d); }
static inline long mpz_matrix_height (mpz_t *A, int d) { return mpz_vec_height (A, d*d); }

static inline long mpz_row_bits (mpz_t *A, int d)
    { long b=0; for ( int i = 0 ; i < d ; i++ ) { b += mpz_height(A[i]); } return b; }
static inline long mpz_matrix_bits (mpz_t *A, int d) { return mpz_row_bits (A, d*d); }

// only reduces if it seems worthwhile
static inline void mpz_vec_mod (mpz_t *A, mpz_t *B, long n, mpz_t m, int *reps)
{
    int i;
    long max, tot, s;

    for ( i = 0,  tot = 0, max = 0 ; i < n ; i++ ) { s = mpz_size (B[i]); if ( s > max ) max = s;  tot += s; }
    if ( 8*max <= 9*mpz_size(m) ) { mpz_vec_set (A, B, n); return ; }   // don't reduce rows whose entries are not much bigger than the modulus
    if ( n == 1 ) { mpz_fdiv_r (A[0], B[0], m); return; }
    s = mpz_mod_fft_crossover(n);
    if ( ! hw_disable_fft && tot > n*s && mpz_size(m) > s ) {
        mpz_vec_mod_fft (A, B, n, m, reps);
    } else {
        mpz_vec_mod_naive (A, B, n, m);
    }
}

static inline void mpz_vec_mod_hard (mpz_t *A, mpz_t *B, long n, mpz_t m)
{
    int i;
    long tot, s;

    for ( i = 0,  tot = 0 ; i < n ; i++ ) { s = mpz_size (B[i]); tot += s; }
    if ( n == 1 ) { mpz_fdiv_r (A[0], B[0], m); return; }
    s = mpz_mod_fft_crossover(n);
    if ( ! hw_disable_fft && tot > n*s && mpz_size(m) > s ) {
        mpz_vec_mod_fft (A, B, n, m, NULL);
    } else {
        mpz_vec_mod_naive (A, B, n, m);
    }
}

static inline void mpz_vec_mod_inplace (mpz_t *A, long n, mpz_t m)
{
    int i;
    long max, tot, s;

    for ( i = 0,  tot = 0, max = 0 ; i < n ; i++ ) { s = mpz_size (A[i]); if ( s > max ) max = s;  tot += s; }
    if ( 8*max <= 9*mpz_size(m) ) return;    // don't reduce rows whose entries are not much bigger than the modulus
    if ( n == 1 ) { mpz_fdiv_r (A[0], A[0], m); return; }
    s = mpz_mod_fft_crossover(n);
    if ( ! hw_disable_fft && tot > n*s && mpz_size(m) > s ) {
        mpz_vec_mod_fft (A, A, n, m, NULL);
    } else {
        mpz_vec_mod_naive (A, A, n, m);
    }
}

static inline void mpz_row_mod (mpz_t *A, mpz_t *B, int d, mpz_t m)
    { mpz_vec_mod (A, B, d, m, NULL); }
    
static inline void mpz_rmatrix_mod (mpz_t *A, mpz_t *B, int r, int d, mpz_t m)
    { mpz_vec_mod (A, B, r*d, m, NULL); }

static inline void mpz_matrix_mod (mpz_t *A, mpz_t *B, int d, mpz_t m)
    { mpz_vec_mod (A, B, d*d, m, NULL); }

static inline void mpz_row_mod_inplace (mpz_t *A, mpz_t *B, int d, mpz_t m)
    { mpz_vec_mod_inplace (A, d, m); }
    
static inline void mpz_matrix_mod_inplace (mpz_t *A, int d, mpz_t m)
    { mpz_vec_mod_inplace (A, d*d, m); }


// computes C = AB, where A is a 1 x d matrix (row vector) and B is a d x d matrix, C cannot alias A!
static inline void mpz_row_matrix_mult (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t w)
{
    if ( d == 1 ) { mpz_mul (C[0], A[0], B[0]); return; }
    long s = d*(long)mpz_row_fft_crossover(d);   // for row-matrix mults require both operands to be larger than crossover
    if ( ! hw_disable_fft && mpz_vec_total_size(A,d)>s && mpz_vec_total_size(B,d*d) > d*s ) {
        mpz_row_matrix_mult_fft (C, A, B, d, w);
    } else {
        mpz_row_matrix_mult_naive (C, A, B, d, w);
    }
}

// computes A = AB, where A is a 1 x d matrix (row vector) and B is a d x d matrix, w should point to d+1 inited mpz_t's
static inline void mpz_row_matrix_mult_inplace (mpz_t *A, mpz_t *B, int d, mpz_t *w)
    { mpz_row_matrix_mult (w, A, B, d, w[d]); mpz_vec_set (A, w, d); }

// computes C = AB, where A is an r x d matrix (r row vectors) and B is a d x d matrix, C cannot alias A or B!
static inline void mpz_rmatrix_mult (mpz_t *C, mpz_t *A, int r, mpz_t *B, int d, mpz_t w, int *reps1, int *reps2)
{
    if ( d == 1 ) { mpz_mul (C[0], A[0], B[0]); return; }
    long s = (r+d)*d*(long)mpz_mat_fft_crossover(d); // we might want to make this depend on r
    if ( ! hw_disable_fft && mpz_vec_total_size(A,r*d) + mpz_vec_total_size(B,d*d) > s ) {
        mpz_rmatrix_mult_fft (C, A, r, B, d, w, reps1, reps2);
    } else {
        mpz_rmatrix_mult_naive (C, A, r, B, d, w);
    }
}

// computes C = AB where A and B are d x d matrices, C cannot alias A or B
static inline void mpz_matrix_mult (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t w, int *reps1, int *reps2) { mpz_rmatrix_mult (C, A, d, B, d, w, reps1, reps2); }

// computes C = AB mod m,  where A is a 1 x d matrix (row vector) and B is a d x d matrix, C cannot alias A!
static inline void mpz_row_matrix_mult_mod (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t m, mpz_t w)
    { mpz_row_matrix_mult (C, A, B, d, w); mpz_row_mod (C, C, d, m); }

// computes A = AB mod m,  where A is a 1 x d matrix (row vector) and B is a d x d matrix, w should point to d+1 inited mpz_t's
static inline void mpz_row_matrix_mult_mod_inplace (mpz_t *A, mpz_t *B, int d, mpz_t m, mpz_t *w)
    { mpz_row_matrix_mult_mod (w, A, B, d, m, w[d]); mpz_vec_set (A, w, d); }

// computes C = AB mod m,  where A is a 1 x d matrix (row vector) and B is a d x d matrix, C cannot alias A!
static inline void mpz_row_matrix_mult_mod_naive (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t m, mpz_t w)
    { mpz_row_matrix_mult_naive (C, A, B, d, w); mpz_vec_mod_naive (C, C, d, m); }

// computes A = AB mod m,  where A is a 1 x d matrix (row vector) and B is a d x d matrix, w should point to d+1 inited mpz_t's
static inline void mpz_row_matrix_mult_mod_naive_inplace (mpz_t *A, mpz_t *B, int d, mpz_t m, mpz_t *w)
    { mpz_row_matrix_mult_mod_naive (w, A, B, d, m, w[d]); mpz_vec_set (A, w, d); }

// computes C = AB mod m, where A and B are d x d matrices, C cannot alias A or B
static inline void mpz_matrix_mult_mod (mpz_t *C, mpz_t *A, mpz_t *B, int d, mpz_t m, mpz_t w)
    { mpz_matrix_mult (C, A, B, d, w, NULL, NULL); mpz_matrix_mod (C, C, d, m); }

// computes C=AB mod m, where A is an r x d matrix (r row vectors) and B is a d x d matrix, C cannot alias A or B!
static inline void mpz_rmatrix_mult_mod (mpz_t *C, mpz_t *A, int r, mpz_t *B, int d, mpz_t m, mpz_t w, int *reps2)
    { mpz_rmatrix_mult (C, A, r, B, d, w, NULL, reps2); mpz_row_mod (C, C, r*d, m); }

// computes A=AB mod m, where A is an r x d matrix (r row vectors) and B is a d x d matrix, w points to r*d+1 mpz_t's
static inline void mpz_rmatrix_mult_mod_inplace (mpz_t *A, int r, mpz_t *B, int d, mpz_t m, mpz_t *w, int *reps2)
    { mpz_rmatrix_mult (w, A, r, B, d, w[r*d], NULL, reps2); mpz_row_mod (A, w, r*d, m); }

// destructively computes the product of  an array of n mpz_t's and puts the result in P
static inline void mpz_vec_product (mpz_t P, mpz_t m[], long n)
{
    if ( n <= 0 ) { mpz_set_one (P); return; }
    if ( n == 1 ) { mpz_set (P, m[0]); return; }
    for ( ; n > 2 ; n = (n+1)/2 ) {
        for ( long i = 0 ; i < n/2 ; i++ ) mpz_mul (m[i],m[i],m[n/2+i]);
        if ( (n&1) ) mpz_set (m[n/2],m[n-1]);
    }
    assert (n==2);
    mpz_mul (P, m[0], m[1]);
}

static inline void mpz_vec_product_ui (mpz_t P, unsigned long m[], long n)
{
    mpz_t *w = mpz_vec_alloc_and_init (n);
    mpz_vec_set_ui(w,m,n); mpz_vec_product (P,w,n);
    mpz_vec_clear_and_free (w, n);
}

static inline mpz_t *mpz_poly_matrix_alloc_and_init (int dim, int deg)
    { return mpz_vec_alloc_and_init (dim*dim*(deg+1)); }
    
static inline void mpz_poly_matrix_clear_and_free (mpz_t *A, int dim, int deg)
    { mpz_vec_clear_and_free (A, dim*dim*(deg+1)); }

// compute g(x) = f(x+c), g cannot overlap f
static inline void mpz_poly_shift_si (mpz_t *g, mpz_t *f, int deg, long c, mpz_t w)
{   
    assert ( f != g );
    assert ( !c || (_si_len(c) + _ui_len(deg)) < 64 );

    mpz_vec_set (g, f, deg+1);
    if ( ! c ) return;
    for ( int d = 1 ; d <= deg ; d++ ) {
        mpz_mul_si (w, f[d], c*d);  //  w = c^1*binom(d,1)
        for ( int i = d-1, j = 2 ; i ; i--, j++ ) {
            mpz_add (g[i], g[i], w); 
            mpz_mul_si (w, w, i*c);  mpz_divexact_ui (w, w, j);    // w = c^j*binomd(d,j)
        }
        mpz_add (g[0], g[0], w);
    }
}

static inline void mpz_poly_shift_up (mpz_t *g, mpz_t *f, int deg, mpz_t w) { mpz_poly_shift_si (g, f, deg, 1, w); }

static inline void mpz_poly_matrix_set (mpz_t *B, mpz_t *A, int dim, int deg) { mpz_vec_set (B, A, dim*dim*(deg+1)); }

// set B[i][j] (x) to A[i][j](x) evaluated at x+1, B cannot overlap A
static inline void mpz_poly_matrix_shift_up (mpz_t *B, mpz_t *A, int dim, int deg, mpz_t w)
    { for ( int i = 0 ; i < dim ; i++ ) for ( int j = 0 ; j < dim ; j++ ) mpz_poly_shift_up (B+(i*dim+j)*(deg+1), A+(i*dim+j)*(deg+1), deg, w); }
static inline void mpz_poly_matrix_shift_si (mpz_t *B, mpz_t *A, int dim, int deg, long c, mpz_t w)
    { for ( int i = 0 ; i < dim ; i++ ) for ( int j = 0 ; j < dim ; j++ ) mpz_poly_shift_si (B+(i*dim+j)*(deg+1), A+(i*dim+j)*(deg+1), deg, c, w); }

// w should point to deg+2 inited mpz_t's
static inline void mpz_poly_matrix_shift_inplace (mpz_t *A, int dim, int deg, long c, mpz_t *w)
{
    for ( int i = 0 ; i < dim ; i++ ) for ( int j = 0 ; j < dim ; j++ )
        { mpz_poly_shift_si (w, A+(i*dim+j)*(deg+1), deg, c, w[deg+1]); mpz_vec_set (A+(i*dim+j)*(deg+1), w, deg+1); }
}

static inline void mpz_poly_mult (mpz_t *h, mpz_t *f, int df, mpz_t *g, int dg)
{
    mpz_vec_set_zero (h, df+dg+1);
    for ( int i = 0 ; i <= df ; i++ ) for ( int j = 0 ; j <= dg ; j++ ) mpz_addmul (h[i+j], f[i], g[j]);
}

// A, B, C are all dim x dim matrices of polys, A and B both degrees degA and degB, while C will have degree degC=degA+degB, w should have space for degC+1 entries
// C cannot overlap A or B
static inline void mpz_poly_matrix_mult (mpz_t *C, mpz_t *A, mpz_t *B, int dim, int degA, int degB, mpz_t *w)
{
    for ( int i = 0 ; i < dim ; i++ ) { // ith row of A
        for ( int j = 0 ; j < dim ; j++ ) {  // jth col of B
            mpz_vec_set_zero (C+(i*dim+j)*(degA+degB+1), degA+degB+1);
            for ( int k = 0 ; k < dim ; k++ ) {  // multiply kth entry in ith row of A by the kth entry in the jth col of B
                mpz_poly_mult (w, A+(i*dim+k)*(degA+1), degA, B+(k*dim+j)*(degB+1), degB);
                mpz_vec_addto (C+(i*dim+j)*(degA+degB+1), w, degA+degB+1);
            }
        }
    }
}

// Compute M(x) = T(x)*T(x+1)*...*T(x+n-1)
// M cannot overlap T
static inline void mpz_poly_matrix_rising_factorial (mpz_t *M, mpz_t *T, int r, int d, int n, mpz_t *w)
{
    int m, m1;

    mpz_poly_matrix_set (M, T, r, d); m = 1;
    if ( m == n ) return;

    // double-and-add
    mpz_t *X = mpz_poly_matrix_alloc_and_init (r, (d*n+1)/2);
    mpz_t *Y = mpz_poly_matrix_alloc_and_init (r, (d*n)/2);
    for ( int i = _ui_len(n)-2 ; i >= 0 ; i-- ) {
        if ( n&(1<<i) ) {
            mpz_poly_matrix_shift_si (Y, T, r, d, m, w[0]); mpz_poly_matrix_mult (X, M, Y, r, m*d, d, w);  m1 = m+1;
        } else {
            mpz_poly_matrix_set (X, M, r, m*d); m1 = m;
        }
        mpz_poly_matrix_shift_si (Y, M, r, m*d, m1, w[0]);
        mpz_poly_matrix_mult (M, X, Y, r, m1*d, m, w);
        m += m1;
    }
    assert (m == n);
    mpz_poly_matrix_clear_and_free (X, r, (d*n+1)/2);
    mpz_poly_matrix_clear_and_free (Y, r, (d*n)/2);
}

// set B[i][j] to A[i][j](x) evaluated at x, given powers of x
static inline void mpz_poly_matrix_eval_pow (mpz_t *B, mpz_t *A, int dim, int deg, mpz_t *xpow)
    { for ( int i = 0 ; i < dim ; i++ ) for ( int j = 0 ; j < dim ; j++ ) mpz_poly_eval_pow (B[i*dim+j], A+(i*dim+j)*(deg+1), deg, xpow); }

// computes Y = M(x), w should point to deg+1 inited mpz_t's
static inline void mpz_poly_matrix_eval_si (mpz_t *Y, mpz_t *M, int dim, int deg, long x, mpz_t *w)
{
    mpz_set_one (w[0]); mpz_set_si (w[1], x);  for ( int i = 2 ; i <= deg ; i++ ) mpz_mul_si (w[i], w[i-1], x);
    for ( int i = 0 ; i < dim ; i++ ) for ( int j = 0 ; j < dim ; j++ ) mpz_poly_eval_pow (Y[i*dim+j], M+(i*dim+j)*(deg+1), deg, w);
}

#ifdef __cplusplus
}
#endif

#endif
