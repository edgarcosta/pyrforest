#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include <math.h>
#include "hwmpz.h"
#include "rtree.h"

// creates a tree with ell+1 levels i=0,...,ell, each of which consists of 2^i vectors of length n (so for d x d matrices, use n = d*d).
mpz_t **rtree_alloc (int ell, int n)
{
    mpz_t **t;
    
    assert ( n >= 1 && ell >= 0 );
    t = hw_malloc ((ell+1)*sizeof(*t));
    for ( int i = 0 ; i <= ell ; i++ ) t[i] = mpz_vec_alloc_and_init (n*(1L<<i));
    return t;
}

// given a tree of d x d matrices with level ell filled in, builds a product tree.
void rtree_build (mpz_t **t, int ell, int d)
{
    mpz_t w;
    assert (ell >= 0 && d > 0);
    if ( d == 1 ) { for ( int i = ell-1 ; i >= 0 ; i-- ) for ( long j = (1L<<i)-1 ; j>=0 ; j-- ) mpz_mul (t[i][j], t[i+1][2*j], t[i+1][2*j+1]);
        return;
    }
    mpz_init2(w,(1<<ell)*mpz_matrix_height(t[ell]+(1<<(ell-1))*d*d,d));
    for ( int i = ell-1 ; i >= 0 ; i-- )  for ( long j = (1L<<i)-1 ; j>=0 ; j-- ) mpz_matrix_mult (t[i]+j*d*d, t[i+1]+2*j*d*d, t[i+1]+(2*j+1)*d*d, d, w);
    mpz_clear(w);
}

// given a tree of d x d matrices with level ell filled in, builds a product tree, removing a known divisor of each product
void rtree_build_div (mpz_t **t, int ell, int d, mpz_t div)
{
    mpz_t w;
    assert (ell >= 0 && d > 0);
    if ( d == 1 ) {
        for ( int i = ell-1 ; i >= 0 ; i-- ) for ( long j = (1L<<i)-1 ; j>=0 ; j-- ) { mpz_mul (t[i][j], t[i+1][2*j], t[i+1][2*j+1]); mpz_divexact (t[i][j], t[i][j], div); }
        return;
    }
    mpz_init2(w,(1<<ell)*mpz_matrix_height(t[ell]+(1<<(ell-1))*d*d,d));
    for ( int i = ell-1 ; i >= 0 ; i-- ) for ( long j = (1L<<i)-1 ; j>=0 ; j-- )
        { mpz_matrix_mult (t[i]+j*d*d, t[i+1]+2*j*d*d, t[i+1]+(2*j+1)*d*d, d, w); mpz_matrix_div (t[i]+j*d*d, d, div); }
    mpz_clear(w);
}

// given a product tree of d x d matrices and a product tree of moduli, and computes the remainder tree of d-vectors starting from V
// w should point to d+1 inited mpz_t's
void rtree_reduce (mpz_t **Rtree, mpz_t *V, mpz_t **Mtree, mpz_t **mtree, int ell, int d, mpz_t *w)
{
    mpz_t *R = w+1;

    assert (ell > 0 && d > 0);

    if ( d == 1 ) { // handle d=1 separately to eliminate overhead of using 1x1 matrices and 1-vectors
        mpz_mod (Rtree[0][0], V[0], mtree[0][0]);
        for (  int i = 1 ; i <= ell ; i++ ) {
            for (  long j = (1L<<i)-1 ; j>=0 ; j-- ) {
                if ( !(j&1) ) {
                    mpz_mod (Rtree[i][j], Rtree[i-1][j/2], mtree[i][j]);
                } else {
                    // don't reduce the sibling matrix (reducing slows things down substantially)
                    mpz_mod (R[0], Rtree[i-1][j/2], mtree[i][j]);
                    mpz_mul (R[0], R[0], Mtree[i][j-1]);
                    mpz_mod (Rtree[i][j], R[0], mtree[i][j]);
                }
            }
        }
    } else {
        int i;
        mpz_row_mod (Rtree[0], V, d, mtree[0][0]);
        for ( i = 1 ; i <= ell-6 ; i++ ) {
            for (  long j = (1L<<i)-1 ; j>=0 ; j-- ) {
                if ( !(j&1) ) {
                    mpz_row_mod (Rtree[i]+j*d, Rtree[i-1]+(j/2)*d, d, mtree[i][j]);
                } else {
                    // don't reduce the sibling matrix (reducing slows things down substantially)
                    mpz_row_mod (R, Rtree[i-1]+(j/2)*d, d, mtree[i][j]);
                    mpz_row_matrix_mult_mod (Rtree[i]+j*d, R, Mtree[i]+(j-1)*d*d, d, mtree[i][j], w[0]);
                }
            }
        }
        // don't bother checking for fft crossovers in the bottom layers of the tree
        for ( ; i <= ell ; i++ ) {
            for (  long j = (1L<<i)-1 ; j>=0 ; j-- ) {
                if ( !(j&1) ) {
                    mpz_vec_mod_naive (Rtree[i]+j*d, Rtree[i-1]+(j/2)*d, d, mtree[i][j]);
                } else {
                    // don't reduce the sibling matrix (reducing slows things down substantially)
                    mpz_vec_mod_naive (R, Rtree[i-1]+(j/2)*d, d, mtree[i][j]);
                    mpz_row_matrix_mult_mod_naive (Rtree[i]+j*d, R, Mtree[i]+(j-1)*d*d, d, mtree[i][j], w[0]);
                }
            }
        }
    }
}

// given a product tree of d x d matrices and a product tree of moduli, and computes the remainder tree of r x d matrices starting from V
void rtree_reduce_rows (mpz_t **Rtree, mpz_t *V, mpz_t **Mtree, mpz_t **mtree, int ell, int d, int r)
{
    mpz_t w;
    mpz_t *R, *M;
     int i;
    
    R = mpz_vec_alloc_and_init (r*d);
    M = mpz_vec_alloc_and_init (d*d);
    mpz_init(w);
    mpz_vec_mod (Rtree[0], V, r*d, mtree[0][0]);
    for ( i = 1 ; i <= ell-4 ; i++ ) {
        for (  long j = (1L<<i)-1 ; j>=0 ; j-- ) {
            if ( !(j&1) ) {
                mpz_vec_mod (Rtree[i]+j*r*d, Rtree[i-1]+(j/2)*r*d, r*d, mtree[i][j]);
            } else {
                mpz_vec_mod (M, Mtree[i]+(j-1)*d*d, d*d, mtree[i][j]);
                mpz_vec_mod (R, Rtree[i-1]+(j/2)*r*d, r*d, mtree[i][j]);
                mpz_rmatrix_mult_mod (Rtree[i]+j*r*d, R, r, M, d, mtree[i][j], w);
            }
        }
    }
    // don't bother checking for fft crossovers in the bottom layers of the tree
    for ( ; i <= ell ; i++ ) {
        for (  long j = (1L<<i)-1 ; j>=0 ; j-- ) {
            if ( !(j&1) ) {
                mpz_vec_mod_naive (Rtree[i]+j*r*d, Rtree[i-1]+(j/2)*r*d, r*d, mtree[i][j]);
            } else {
                mpz_vec_mod_naive (R, Rtree[i-1]+(j/2)*r*d, r*d, mtree[i][j]);
                // don't reduce the sibling matrix
                mpz_rmatrix_mult_mod (Rtree[i]+j*r*d, R, r, Mtree[i]+(j-1)*d*d, d, mtree[i][j], w);
            }
        }
    }
    mpz_clear (w);
    mpz_vec_clear_and_free (R, r*d);
    mpz_vec_clear_and_free (M, d*d);
}

void rtree_free (mpz_t **t, int ell, int n)
    {  for ( int i = 0 ; i <= ell ; i++ ) { mpz_vec_clear_and_free (t[i], n*(1L<<i)); } hw_free (t, (ell+1)*sizeof(*t)); }
