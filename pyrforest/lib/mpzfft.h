/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_MPZFFT_H
#define ZZ_MPZFFT_H

#include "mpnfft.h"
#include "mpnfft_mod.h"


typedef zz_mpnfft_params_t mpzfft_params_t;
typedef zz_mpnfft_poly_t mpzfft_t;


// initialise params for use with given bits, terms, num_primes, moduli
// moduli must stay in scope for the life of the object
static inline
void mpzfft_params_init(mpzfft_params_t* params, size_t bits, unsigned terms,
			unsigned num_primes, zz_moduli_t* moduli)
{
  zz_mpnfft_params_init(params, bits, terms, num_primes, moduli);
}

// destroy params
static inline void mpzfft_params_clear(mpzfft_params_t* params)
{
  zz_mpnfft_params_clear(params);
}

// initialise op for use with params
// params must stay in scope for the life of op
// value is set to zero
// no memory is allocated yet; that happens as needed
static inline void mpzfft_init(mpzfft_t op, mpzfft_params_t* params)
{
  zz_mpnfft_poly_init(op, params);
}

// destroy op
static inline void mpzfft_clear(mpzfft_t op)
{
  zz_mpnfft_poly_clear(op);
}

// free up space occupied by op (without destroying it), value set to zero
static inline void mpzfft_dealloc(mpzfft_t op)
{
  zz_mpnfft_poly_dealloc(op);
}


// convert op to FFT representation in rop
void mpzfft_fft(mpzfft_t rop, mpz_t op, int threads);

// convert op out of FFT representation, write result to rop
// NOTE: the contents of op are DESTROYED, and op is set to zero
void mpzfft_ifft(mpz_t rop, mpzfft_t op, int threads);

// rop := op
// rop, op must have same params
// parameter aliasing allowed
static inline void mpzfft_set(mpzfft_t rop, mpzfft_t op, int threads)
{
  zz_mpnfft_poly_set(rop, op, threads);
}

// rop := op1 + op2
// rop, op1, op2 must have same params
// parameter aliasing allowed
static inline void mpzfft_add(mpzfft_t rop, mpzfft_t op1, mpzfft_t op2,
			      int threads)
{
  zz_mpnfft_poly_add(rop, op1, op2, threads);
}


// rop := op1 - op2
// rop, op1, op2 must have same params
// parameter aliasing allowed
static inline void mpzfft_sub(mpzfft_t rop, mpzfft_t op1, mpzfft_t op2,
			      int threads)
{
  zz_mpnfft_poly_sub(rop, op1, op2, threads);
}


// rop := -op1
// rop, op must have same params
// parameter aliasing allowed
static inline void mpzfft_neg(mpzfft_t rop, mpzfft_t op, int threads)
{
  zz_mpnfft_poly_neg(rop, op, threads);
}


// rop := op1 * op2
// rop, op1, op2 must have same params
// parameter aliasing allowed
static inline void mpzfft_mul(mpzfft_t rop, mpzfft_t op1, mpzfft_t op2,
			      int threads)
{
  zz_mpnfft_poly_mul(rop, op1, op2, 0, threads);
}

// multiplies matrices of fourier coefficients, rop := op1 * op2
// op1: dim1 rows, dim2 columns
// op2: dim2 rows, dim3 columns
// rop: dim1 rows, dim3 columns
// expects that params has terms >= dim2
// matrices are stored in row-major order
// parameter aliasing allowed
static inline
void mpzfft_matrix_mul(mpzfft_t* rop, mpzfft_t* op1, mpzfft_t* op2,
		       unsigned dim1, unsigned dim2, unsigned dim3,
		       int threads)
{
  zz_mpnfft_poly_matrix_mul(rop, op1, op2, dim1, dim2, dim3, threads);
}



typedef zz_mpnfft_mod_t mpzfft_mod_t;


// initialise "mod" for computing x mod d, where |x| has at most "bits" bits,
// d > 0. Must have bits >= number of bits in d.
void mpzfft_mod_init(mpzfft_mod_t* mod, size_t bits, mpz_t d,
		     unsigned num_primes, zz_moduli_t* moduli, int threads);


static inline void mpzfft_mod_clear(mpzfft_mod_t* mod)
{
  zz_mpnfft_mod_clear(mod);
}


// rop := op mod d, with 0 <= rop < d
// |op| must be at most n bits (any sign allowed)
// rop may alias op
void mpzfft_mod_mod(mpzfft_mod_t* mod, mpz_t rop, mpz_t op, int threads);


#endif
