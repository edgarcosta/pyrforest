#ifndef _RTREE_INCLUDE_
#define _RTREE_INCLUDE_

#include <stdio.h>

#define RTREE_MAX_LEVELS    30

mpz_t **rtree_alloc (int ell, int n);   // n = # entries in each element, either r*r or rows*r
void rtree_free (mpz_t **t, int ell, int n);

/*
  Input:  M = rtree_alloc_init(ell,r*r) with M[ell][j] := M_j filled with r x r inputs for 0 <= j < 2^ell
  Output: M[i,j] = M[i+1][2*j] * M[i+1,j*2+1] for 0 <= i < ell and 0 <= j < 2^i
*/
void rtree_build (mpz_t **Mtree, int ell, int r);
void rtree_build_div (mpz_t **Mtree, int ell, int r, mpz_t div);    // div is an integer that will be removed from each pair-wise product (must exactly divide)

/*
  Input:  M = rtree_alloc(ell,r*r);  rtree_build(M,ell,r), m = rtree_alloc(r,1);  rtree_build(m,ell,r);  R = rtree(alloc(ell,r);
  Output: R[0][0] = V mod m[0][0]
          R[i,j] = R[i-1][j/2] mod m[i][j] (even j), R[i-1][j/2] = R[i-1][j/2]*M[i,j-1] mod m[i][j] (odd j), for 1 <= i <= ell, 0 <= j < 2^i
          In particular, R_j := R[ell,j] = V*M_0*M_1...*M_{j-1} mod m_j := m[ell,j]
  w should point to r+1 inited mpz_t's
*/
void rtree_reduce (mpz_t **Rtree, mpz_t *V, mpz_t **Mtree, mpz_t **mtree, int ell, int r, mpz_t *w);

// Same as rtree_reduce except now V is a rows*r matrix instead of a single row vector of length r
void rtree_reduce_rows (mpz_t **Rtree, mpz_t *V, mpz_t **Mtree, mpz_t **mtree, int ell, int r, int rows, int *reps);

void rtree_remainders_d (mpz_t V[], mpz_t z, mpz_t C[], mpz_t M[], int d, mpz_t **mtree, int ell);

#endif
