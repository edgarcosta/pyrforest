#include <assert.h>
#include "hwmem.h"
#include "hwmpz.h"
#include "rtree.h"
#include "rforest.h"
#include "time.h"

#define _max(a,b) ((a)>(b)?(a):(b))

void rforest (mpz_t *A, mpz_t *V, int rows, mpz_t *M, int deg, int dim, mpz_t *m, long kbase, long *ks, long n, mpz_t z, int kappa, int *reps)
{
    assert ( A && V && rows > 0 && M && deg >= 0 && dim > 0 && m && ks && n >= 0 && kappa >= 0 );
    if ( !n ) return;

    int ell = _ui_len(n) - kappa; // we will use <= 2^kappa trees of height ell.
    if ( ell < 0 ) ell = 0;

    hw_mem_init(0);
    hw_mpz_setup();

    // working space
    mpz_t *w = mpz_vec_alloc_and_init2 (_max(deg,rows*dim)+1, mpz_bits(z));

    if ( !ell ) {

        mpz_t *Mk = mpz_vec_alloc_and_init (dim*dim);
        for ( long k = kbase, i = 0 ; i < n ; i++ ) {
            long nextk = ks[i];  assert(nextk >= k);
            while ( k < nextk ) {
                mpz_poly_matrix_eval_si (Mk, M, dim, deg, k++, w);
                mpz_rmatrix_mult_mod_inplace (V, rows, Mk, dim, z, w, 0);
            }
            mpz_vec_mod_hard (A+i*rows*dim, V, rows*dim, m[i]);
            mpz_divexact (z, z, m[i]);
            mpz_rmatrix_mod (V, V, rows, dim, z);
        }
        mpz_vec_clear_and_free (Mk, dim*dim);

    } else {

        mpz_t *W0 = mpz_matrix_alloc_and_init (dim), *W1 = mpz_matrix_alloc_and_init (dim);
        mpz_t **mtree = rtree_alloc (ell, 1), **Mtree = rtree_alloc (ell, dim*dim), **Rtree = rtree_alloc (ell, rows*dim);
        long i, s, k = kbase, t = 1L << ell;
        clock_t seconds = clock();
        for ( s = 0 ; s <= n ; s += t ) {   // we need to insert 1 at the start of the modulus array so we have n+1 moduli
            long i0 = s;
            long i1 = i0+t;

            // set moduli for current subtree, pad with 1's if necessary, note that i shifts by 1
            if ( !i0 ) mpz_set_one(mtree[ell][0]);  // insert 1 at the begining
            for ( long i = i0?i0:1 ; i < i1 ; i++ ) if ( i <= n ) mpz_set (mtree[ell][i-i0],m[i-1]); else mpz_set_one(mtree[ell][i-i0]);

            // set leaves of matrix tree using transition matrix
            // note that i is not shifted by 1
            mpz_t *Mk = Mtree[ell];
            for ( i = i0 ; i < i1 && i < n ; i++, Mk += dim*dim ) {
                long nextk = ks[i];  assert(nextk >= k);
                if ( k == nextk ) mpz_matrix_set_one (Mk, dim); else mpz_poly_matrix_eval_si (Mk, M, dim, deg, k++, w);
                while ( k < nextk ) {
                    mpz_matrix_set (W0, Mk, dim);
                    mpz_poly_matrix_eval_si (W1, M, dim, deg, k++, w);
                    mpz_matrix_mult (Mk, W0, W1, dim, w[0]);
                }
            }
            for ( ; i < i1 ; i++, Mk += dim*dim ) mpz_matrix_set_one (Mk, dim); // pad with the identity matrix

            printf("a %d\n", clock()-seconds);
            seconds = clock();
            // Run the remainder tree algorithm: build/build/reduce
            rtree_build (mtree, ell, 1);  rtree_build (Mtree, ell, dim);
            printf("aa %d\n", clock()-seconds);
            seconds = clock();
            rtree_reduce_rows (Rtree, V, Mtree, mtree, ell, dim, rows);
            printf("b %d\n", clock()-seconds);
            seconds = clock();

            // copy leaves of current tree to output -- use mpz_vec_mod_naive to force a hard mod (leaf values should be small)
            // note that i shifts by 1
            for ( i = i0?i0:1 ; i < i1 && i <= n ; i++ ) {
                mpz_vec_mod_naive (A+(i-1)*rows*dim, Rtree[ell] + (i-i0)*(rows*dim), rows*dim, m[i-1]);
            }
            // update the modulus and transfer vector
            mpz_divexact (z, z, mtree[0][0]);
            printf("c %d\n", clock()-seconds);
            seconds = clock();
            mpz_rmatrix_mult_mod_inplace (V, rows, Mtree[0], dim, z, w, reps);
            printf("d %d\n", clock()-seconds);
            seconds = clock();
        }
        mpz_matrix_clear_and_free (W0, dim); mpz_matrix_clear_and_free (W1, dim);
        rtree_free (mtree, ell, 1);  rtree_free (Mtree, ell, dim*dim);  rtree_free (Rtree, ell, rows*dim);
    }

    mpz_vec_clear_and_free (w, _max(deg,rows*dim)+1);
    hw_mpz_clear();
    hw_mem_clear();

    // for consistency reduce V = V mod z before returning (note mpz_rmatrix_mod does a soft reduction)
    mpz_vec_mod_naive (V, V, rows*dim, z);
}


void mproduct (mpz_t z, mpz_t *m, long n)
{
    mpz_t *w = mpz_vec_alloc_and_init (n);
    mpz_vec_set (w, m, n);
    mpz_vec_product (z, w, n);
}
