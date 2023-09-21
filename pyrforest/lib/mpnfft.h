/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_MPNFFT_H
#define ZZ_MPNFFT_H

#include "mpzfft_moduli.h"


typedef struct
{
  zz_moduli_t* moduli;
  unsigned num_primes;
  size_t terms;
  unsigned r;

  // FFT transform length = N = 2^lgN
  unsigned lgN;
  size_t N;

  // points = number of fourier transform points used; this is guaranteed to be
  // enough to recover integers whose number of bits is "bits" as specified in
  // the call to zz_mpnfft_params_init()
  // (must be an admissible transform size <= N, see fft62.h)
  size_t points;
}
zz_mpnfft_params_t;


// chooses parameters so that the transform can accommodate any sum of products
// such that |x_1*y_1| + ... + |x_{terms}*y_{terms}| < 2^bits.
void zz_mpnfft_params_init(zz_mpnfft_params_t* params, size_t bits,
			   size_t terms, unsigned num_primes,
			   zz_moduli_t* moduli);

static inline void zz_mpnfft_params_clear(zz_mpnfft_params_t* params) { ; }


/*
  Represents a polynomial in Z[x].

  There are two possible representations: coefficient representation and
  FFT representation.

  In coefficient representation, coefficients are represented by their residues
  moduli several FFT primes q_i.

  In FFT representation, each f(x) modulo q_i is represented by its values at
  exactly "points" roots of unity in Z/q_i.
*/
typedef struct
{
  zz_mpnfft_params_t* params;

  // In coefficient representation, size = number of meaningful coefficients in
  // each data[i]; remaining coefficients are assumed zero.
  // Must always have 0 <= size <= points.
  // If size is nonzero, it must be an admissible size (see fft62.h),
  // In FFT representation, either size == 0 (represents zero), or
  // size == points.
  size_t size;

  // If data[0] == NULL, then no space is allocated yet.
  // Otherwise, data[i] = buffer of length "points" = coefficients for i-th
  // prime, in [0, 2p), for 0 <= i < num_primes.
  uint64_t* data[ZZ_MAX_PRIMES];
}
zz_mpnfft_poly_struct;

typedef zz_mpnfft_poly_struct zz_mpnfft_poly_t[1];


// initial value is zero
void zz_mpnfft_poly_init(zz_mpnfft_poly_t P, zz_mpnfft_params_t* params);

void zz_mpnfft_poly_clear(zz_mpnfft_poly_t P);

// allocate memory if not already allocated
void zz_mpnfft_poly_alloc(zz_mpnfft_poly_t P);

// frees memory without destroying object
void zz_mpnfft_poly_dealloc(zz_mpnfft_poly_t P);

// Converts integer {up,un} to polynomial. Any un >= 0 is ok.
// result is multiplied by sign * 2^lgS, where sign = 1 or -1.
// additionally if scale == 1, result is multiplied by 1/N mod p.
// (note that all this scaling happens essentially for free because it gets
// folded into the split_reduce step)
// Writes as many coefficients as necessary (always an admissible size or
// zero), up to a maximum of "points".
void zz_mpnfft_mpn_to_poly(zz_mpnfft_poly_t P, mp_limb_t* up, size_t un,
			   int sign, int lgS, int scale, int threads);

// Writes f(2^r) mod B^rn to {rp,rn}
void zz_mpnfft_poly_to_mpn(mp_limb_t* rp, size_t rn, zz_mpnfft_poly_t P,
			   int threads);

// rop := op
// (works for both coefficient and FFT representation)
void zz_mpnfft_poly_set(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op, int threads);

// rop := -op
// (works for both coefficient and FFT representation)
void zz_mpnfft_poly_neg(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op, int threads);

// rop := op1 + op2
// (works for both coefficient and FFT representation)
void zz_mpnfft_poly_add(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op1,
			zz_mpnfft_poly_t op2, int threads);

// rop := op1 - op2
// (works for both coefficient and FFT representation)
void zz_mpnfft_poly_sub(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op1,
			zz_mpnfft_poly_t op2, int threads);

// computes FFT of op (zero or "points" points), writes result to rop
// rop may alias op
void zz_mpnfft_poly_fft(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op,
			int threads);

// computes IFFT of op, writes result to rop (zero or "points" coefficients)
// rop may alias op
// if scale == 1, result is multiplied by 1/N mod p
void zz_mpnfft_poly_ifft(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op,
			 int scale, int threads);

// multiply fourier coefficients pointwise, rop := op1 * op2
// all parameters may alias each other
// if use_pinvb == 1, it uses pinvb (montgomery) multiplication, so the
// result ends up multiplied by 1/2^64 mod p.
void zz_mpnfft_poly_mul(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op1,
			zz_mpnfft_poly_t op2, int use_pinvb, int threads);

// multiplies matrices of fourier coefficients, rop := op1 * op2
// op1: dim1 rows, dim2 columns
// op2: dim2 rows, dim3 columns
// rop: dim1 rows, dim3 columns
// expects that params has "terms" >= dim2
// matrices are stored in row-major order
// parameter aliasing allowed
void zz_mpnfft_poly_matrix_mul(zz_mpnfft_poly_t* rop,
			       zz_mpnfft_poly_t* op1,
			       zz_mpnfft_poly_t* op2,
			       unsigned dim1, unsigned dim2, unsigned dim3,
			       int threads);


#endif
