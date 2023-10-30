#ifndef _INCLUDE_RFOREST_
#define _INCLUDE_RFOREST_

#ifdef __cplusplus
extern "C" {
#endif


/*
    Given a transiion matrix over Z[x]^(dim x dim) with entries
    of degree at most deg, a list of n moduli m, and an initial
    row x dim matrix V, rforest computes, for 0 <= i < nL

        A[i] = V * prod_j=kbase^(k[i]-1) M(j) mod m[i]

    where A[i] denotes a row x dim matrix stored at offset
    i*(row*dim) in the array A, and sets

        V = V*prod_j=kbase^k[n-1]-1)M(j)
        z = z / prod_i=0^(n-1) m[i]
*/
void rforest (mpz_t *A,   // array of size rows*dim*n (outputs)
              mpz_t *V,   // array of rows*dims (updated on return)
              int rows,   // number of rows of V
              mpz_t *M,   // array of size dim*dim*(deg+1)
              int deg,    // max deg of polynomial entries of M
              int dim,    // dimension of square matrix M over Z[x]
              mpz_t *m,   // array of n positive moduli
              long kbase, // first value of k
              long *k,    // array of n values of k range endpoints
              long n,     // number of moduli and rows*dim outputs
              mpz_t z,    // integer divisible by product of the moduli
              int kappa,  // log_2 of number of trees in the forest
              int *reps); // specify repeated entries of products (or NULL)

// computes z = prod_i=0^(n-1) m[i] using a product tree
void mproduct (mpz_t z, mpz_t *m, long n);

#ifdef __cplusplus
}
#endif

#endif
