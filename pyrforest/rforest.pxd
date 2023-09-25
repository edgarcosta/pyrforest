from sage.libs.gmp.types cimport mpz_t

cdef extern from "<rforest.h>":
    void mproduct (mpz_t z, mpz_t *m, long n)

    void rforest (mpz_t *A,   # array of size rows*dim*n (outputs)
                  mpz_t *V,   # array of rows*dim (updated on return)
                  int rows,   # number of rows of V
                  mpz_t *M,   # array of size dim*dim*(deg+1)
                  int deg,    # max deg of polynomial entries of M
                  int dim,    # dimension of square matrix M over Z[x]
                  mpz_t *m,   # array of n positive moduli
                  long kbase, # initial offset
                  long *k,    # array of n values of k
                  long n,     # number of moduli and rows*dim outputs
                  mpz_t z,    # integer divisible by product of the moduli
                  int kappa)  # log_2 of number of trees in the forest
