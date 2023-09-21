/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_FERMAT_H
#define ZZ_FERMAT_H

#include "mpzfft_moduli.h"


/*
  A fermat modulus is a product of the form
      F = (B^{t*M_0} + 1) * ... * (B^{t*M_{w-1}} + 1)
  where M_0 > ... > M_{w-1} are distinct powers of two, and t >= 1. Here
  w = Hamming weight of M = M_0 + ... + M_{w-1}. These parameters are stored
  in a zz_fermat_params_t object.

  Let u be in [0, F). The fermat transform of u is defined to be the vector of
  residues
     u_i  mod  B^{t*M_i} + 1,     for 0 <= i < w.
  Such a vector is represented by a zz_fermat_transform_t object. The functions
     zz_fermat_mpn_to_transform
     zz_fermat_mpz_to_transform
     zz_fermat_transform_to_mpn
     zz_fermat_transform_to_mpz
  convert to and from the fermat transform.

  A fermat polynomial is a vector of integer polynomials
     U_i(x)  mod  x^{K*M_i} + 1,    for 0 <= i < w,
  where K = 2^lgK is a power of two. Such a vector is represented by a
  zz_fermat_poly_t object. The (signed) coefficients are stored modulo several
  FFT primes; the number of primes used implies a bound on the absolute value of
  the coefficients.

  The function
     zz_fermat_transform_to_poly
  converts a fermat transform to a fermat polynomial, by splitting each u_i into
  chunks of r bits. Thus we get
     U_i(2^r) = u_i  mod  B^{t*M_i} + 1.
  We require r*K = 64*t, so that the evaluation map makes sense. The function
     zz_fermat_poly_to_transform
  goes the other way, performing the evaluation
     U_i -> U_i(2^r) = u_i  mod  B^{t*M_i} + 1.

  The zz_fermat_poly_t object can also represent U_i(x) in FFT representation,
  i.e. by storing the values of the U_i(x) at the roots of each x^{K*M_i} + 1
  modulo each FFT prime. The corresponding evaluation/interpolation maps are
  performed by the functions
     zz_fermat_poly_fft
     zz_fermat_poly_ifft

  The point of all this is of course to be able to perform efficient
  multiplication modulo F. The parameters are chosen so that the bound for the
  coefficients for the U_i permits sums of up to "terms" products U * V before
  converting out of FFT representation.
*/

#define ZZ_FERMAT_MAX_WEIGHT 8
#define ZZ_FERMAT_PADDING (ZZ_MAX_PRIMES + 3)

#if TEST
#define ZZ_FERMAT_SLICE 4
#else
#define ZZ_FERMAT_SLICE 1024
#endif


typedef struct
{
  zz_moduli_t* moduli;
  unsigned num_primes;
  size_t terms;
  size_t M;
  unsigned w;
  unsigned e[ZZ_FERMAT_MAX_WEIGHT];   // M_i = 2^e[i], 0 <= i < w
  unsigned r;
  size_t t;
  unsigned lgK;
  size_t K;
}
zz_fermat_params_t;


// select fermat modulus such that F >= B^n.
void zz_fermat_params_init(zz_fermat_params_t* params, size_t n, size_t terms,
			   unsigned num_primes, zz_moduli_t* moduli);

static inline void zz_fermat_params_clear(zz_fermat_params_t* params) { ; }

// number of limbs needed for the main buffer of zz_fermat_transform_t
static inline size_t zz_fermat_params_bufsize(zz_fermat_params_t* params)
{
  return params->M * params->t + ZZ_FERMAT_PADDING;
}


// writes fermat modulus F to {rp, t*M+1}
void zz_fermat_params_F_mpn(zz_fermat_params_t* params, mp_limb_t* rp);

// writes fermat modulus F to rop
void zz_fermat_params_F_mpz(zz_fermat_params_t* params, mpz_t rop);


/*
  Stores a vector of fermat residues in the following format.

  The main buffer ("data") has length t * M. It contains the low t * M_0 limbs
  of u_0, followed by the low t * M_1 limbs of u_1, ..., and finally the low
  t * M_{w-1} limbs of u_{w-1}. The high limb of each u_i is stored in hi[i],
  which is a *signed* value in {-1, 0, 1}, i.e. the u_i are not maintained in
  the canonical interval [0, B^(t*M_i) + 1); instead they are allowed to float
  around in the larger interval [-B^(t*M_i), 2*B^(t*M_i)).
*/
typedef struct
{
  zz_fermat_params_t* params;

  // data = main buffer of length t * M + ZZ_FERMAT_PADDING.
  // if own == 1, then buffer was allocated by zz_fermat_tranform_init() and
  //    will be freed by zz_fermat_transform_clear().
  // if own == 0, then buffer was provided by caller of
  //    zz_fermat_transform_init(), and must be freed by the caller.
  int own;
  mp_limb_t* data;

  mp_limb_signed_t hi[ZZ_FERMAT_MAX_WEIGHT];
}
zz_fermat_transform_struct;

typedef zz_fermat_transform_struct zz_fermat_transform_t[1];


/* Initialises T with given parameters.
   If data == NULL, buffer is allocated internally.
   Otherwise data must of size returned by zz_fermat_params_bufsize().
   Initial value of residues is undefined. */
void zz_fermat_transform_init(zz_fermat_transform_t T,
			      zz_fermat_params_t* params, mp_limb_t* data);

void zz_fermat_transform_clear(zz_fermat_transform_t T);


// adds c*F to {op, t*M+1}
void zz_fermat_add_F(mp_limb_t* op, mp_limb_signed_t c,
		     zz_fermat_params_t* params);


/* Computes fermat transform of {op,n}. Any any >= 0 is ok.
   Inplace operation is permitted, i.e. if T was initialised with T->data == op,
   provided of course that there is enough room in op. */
void zz_fermat_mpn_to_transform(zz_fermat_transform_t T,
				mp_limb_t* op, size_t n, int threads);

// Computes fermat transform of op (op may be positive, negative or zero).
void zz_fermat_mpz_to_transform(zz_fermat_transform_t T, mpz_t op, int threads);

/* Combines fermat residues to obtain integer in [0, F).
   Writes exactly t * M + 1 limbs.
   Inplace operation is permitted. */
void zz_fermat_transform_to_mpn(mp_limb_t* rop, zz_fermat_transform_t T,
				int threads);

// Combines fermat residues to obtain integer in [0, F).
void zz_fermat_transform_to_mpz(mpz_t rop, zz_fermat_transform_t T,
				int threads);


// rop := 0
void zz_fermat_transform_zero(zz_fermat_transform_t T, int threads);


// rop := op1 - op2
// aliasing allowed
void
zz_fermat_transform_sub(zz_fermat_transform_t rop, zz_fermat_transform_t op1,
			zz_fermat_transform_t op2, int threads);

// rop := -op
// aliasing allowed
void
zz_fermat_transform_neg(zz_fermat_transform_t rop, zz_fermat_transform_t op,
			int threads);



typedef struct
{
  zz_fermat_params_t* params;

  // data[i][j] = data for j-th prime for the component modulo x^{M_i*K} + 1.
  // Coefficients in [0, 2p).
  // If data[0][0] == NULL, then no space is allocated yet.
  // Otherwise, data[i][j] has length M_i * K, for 0 <= i < w,
  // 0 <= j < num_primes.
  uint64_t* data[ZZ_FERMAT_MAX_WEIGHT][ZZ_MAX_PRIMES];
}
zz_fermat_poly_struct;

typedef zz_fermat_poly_struct zz_fermat_poly_t[1];


// initial value of polynomial is undefined
void zz_fermat_poly_init(zz_fermat_poly_t P, zz_fermat_params_t* params);

void zz_fermat_poly_clear(zz_fermat_poly_t P);

// ensures that space is allocated for polynomials
void zz_fermat_poly_alloc(zz_fermat_poly_t P);

// frees memory without destroying object
void zz_fermat_poly_dealloc(zz_fermat_poly_t P);


// convert from fermat transform to fermat polynomial
// result is multiplied by sign * 2^lgS, where sign = 1 or -1.
// additionally if scale == 1, the component modulo x^{M_i*K} + 1 is multiplied
// by 1/(M_i*K) mod p.
// (note that all this scaling happens essentially for free because it gets
// folded into the split_reduce step)
void zz_fermat_transform_to_poly(zz_fermat_poly_t P, zz_fermat_transform_t T,
				 int sign, int lgS, int scale, int threads);

void zz_fermat_poly_to_transform(zz_fermat_transform_t T, zz_fermat_poly_t P,
				 int threads);

// computes FFT of op, writes result to rop
// rop may alias op
void zz_fermat_poly_fft(zz_fermat_poly_t rop, zz_fermat_poly_t op,
			int threads);

// computes IFFT of op, writes result to rop
// rop may alias op
// note: this performs no scaling, i.e. component modulo x^{M_i*K} + 1
// is multiplied by M_i*K mod p; this can be undone by using the scaling
// parameter in zz_fermat_transform_to_poly()
void zz_fermat_poly_ifft(zz_fermat_poly_t rop, zz_fermat_poly_t op,
			 int threads);

// rop := op1 + op2, adds coefficients pointwise
void zz_fermat_poly_add(zz_fermat_poly_t rop, zz_fermat_poly_t op1,
			zz_fermat_poly_t op2, int threads);

// rop := op1 - op2, subtracts coefficients pointwise
void zz_fermat_poly_sub(zz_fermat_poly_t rop, zz_fermat_poly_t op1,
			zz_fermat_poly_t op2, int threads);

// multiply fourier coefficients pointwise, rop := op1 * op2
// if use_pinvb == 1, it uses pinvb (montgomery) multiplication, so the
// result ends up multiplied by 1/2^64 mod p.
void zz_fermat_poly_mul(zz_fermat_poly_t rop, zz_fermat_poly_t op1,
			zz_fermat_poly_t op2, int use_pinvb, int threads);


#endif
