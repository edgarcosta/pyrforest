#include <assert.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "rforest.h"

static inline double get_time (void)
    { struct timeval time;  gettimeofday(&time, NULL);  return time.tv_sec + time.tv_usec / 1000000.0; }

static long aplift (long ap, long p)
    { return ap > p/2 ? ap-p : ( ap < -p/2 ? ap+p : ap); }

int main (int argc, char *argv[])
{
    mpz_t *A, *V, *M, *m;
    mpz_t z, w, D;
    long a, b, i, n, N, *k;
    int rows, dim, deg, kappa;

    if ( argc < 5 ) {
        puts ("usage: test_rforest a b N kappa");
        puts ("       Computes sum of Frobenius traces a_p of y^2 = x^3 + ax + b for p in [17,N] not dividing 4A^3+27B^2\n");
        return 0;
    }
    a = atol(argv[1]); b = atol(argv[2]); N = atol(argv[3]); kappa = atoi(argv[4]);
    if ( ! b ) { printf ("Please make b nonzero.\n"); return -1; }
    if ( N < 17 ) { printf ("N=%ld should be at least 17\n", N); return -1; }
    mpz_init (D); mpz_init (w);
    mpz_set_si (w, b); mpz_mul_si (w, w, b);  mpz_mul_ui (w, w, 27);
    mpz_set_si (D, a); mpz_mul_si (D, D, a); mpz_mul_si (D, D, a);  mpz_mul_2exp (D, D, 2);  mpz_add (D, D, w);
    if ( ! mpz_sgn(D) ) { printf ("y^2 = x^3 %+ld*x %+ld is singular\n", a, b); return -1; }

    // stupidly enumerate primes
    printf ("Stupidly enumerating primes in [17,%ld]...", N);
    long maxn = ceil (1.25506*N / log(N)) + 1; // upper bound on pi(N)
    m = malloc (maxn * sizeof(*m));
    mpz_init_set_ui(m[0], 17);
    for ( n = 1 ; n < maxn; n++ ) {
        mpz_init (m[n]); mpz_nextprime (m[n], m[n-1]);
        if ( mpz_cmp_ui (m[n], N) > 0 ) { mpz_clear(m[n]); break; }
    }
    assert (n < maxn);
    mpz_init2 (z ,N);
    mproduct (z, m, n);
    printf ("%ld primes, product has %ld bits.\n", n, mpz_sizeinbase(z,2));

    rows = deg = 1; dim = 3;
    V = malloc (sizeof(*V)*rows*dim); for ( i = 0 ; i < rows*dim ; i++ ) mpz_init (V[i]);
    mpz_set_ui (V[dim-1], 1);
    A = malloc (sizeof(*A)*rows*dim*n); for ( i = 0 ; i < n*rows*dim ; i++ ) mpz_init (A[i]);

    // Create transition matrix for y^2=cubic with nonzero constant term
    M = malloc ((deg+1)*dim*dim*sizeof(*M));
    for ( i = 0 ; i < (deg+1)*dim*dim ; i++ ) mpz_init (M[i]);
    mpz_set_ui (M[4], 3);
    mpz_set_si(M[5], -2);        // M[0,2] = (3-2*k)*f_3 = 3 - 2*k
    mpz_set_si (M[7], 2*b);      // M[1,0] = 2*k*f_0 = 2*b*k, note M[1,2] = (2-2k)f_2 = 0 
    mpz_set_si (M[15], 2*b);     // M[2,1] = 2*k*f+0 = 2*b*k
    mpz_set_si (M[16], a);
    mpz_set_si (M[17], -2*a);    // M[2,2] = (1-2*k)*f_1 = a - 2*a*k
 
    k = malloc(n*sizeof(*k));
    for ( i = 0 ; i < n ; i++ ) k[i] = mpz_get_ui(m[i]); // for m[i]=p we want to the product of M(k) over k < p

    printf ("Calling rforest with rows=%d, deg=%d, dim=%d, n=%ld, kappa=%d\n", rows, deg, dim, n, kappa);
    double t = get_time();
    rforest (A, V, rows, M, deg, dim, m, 1, k, n, z, kappa);
    printf ("Took %.3fs\n", get_time()-t);

    // sanity check that A values are all reduced
    for ( i = 0 ; i < n ; i++ ) for ( int j = 0 ; j < rows*dim ; j++ ) assert (mpz_cmp(A[i*rows*dim+j],m[i]) < 0);

    // sanity check that z=1 and V=0 on return
    assert (mpz_cmp_ui(z,1) == 0);
    for ( i = 0 ; i < rows*dim ; i++ ) assert (mpz_sgn(V[i]) == 0);

    // compute sum of good a_p for p in [17,N] by lifting last entry in each A_i to Z
    long sum = 0;
    mpz_set_si(w,b);
    for ( i = 0 ; i < n ; i++ ) {
        long p = mpz_get_ui(m[i]);
        if ( ! mpz_divisible_ui_p (D, p) ) {
            long ap = aplift (- mpz_kronecker_ui(w,p) * mpz_get_ui(A[3*i+2]), p);
            sum += ap;
        }
    }

    printf ("Trace sum over good p in [17,%ld] is %ld\n", N, sum);

    free(k);
    for ( i = 0 ; i < (deg+1)*dim*dim ; i++ ) mpz_clear(M[i]);
    for ( i = 0 ; i < rows*dim ; i++ ) mpz_clear(V[i]);
    for ( i = 0 ; i < n*rows*dim ; i++ ) mpz_clear(A[i]);
    for ( i = 0 ; i < n ; i++ ) mpz_clear (m[i]);
    mpz_clear (w); mpz_clear(z);
    free (m);
    return 0;
}
