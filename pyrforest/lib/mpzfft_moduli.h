/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_MODULI_H
#define ZZ_MODULI_H

#include <gmp.h>
#include "fft62/fft62.h"


// maximum number of FFT primes
#define ZZ_MAX_PRIMES 8


// global table of usable FFT primes
extern uint64_t global_p[ZZ_MAX_PRIMES];


/*
  Precomputed data for a list of FFT primes
*/
typedef struct
{
  // number of primes, must have 1 <= num_primes <= ZZ_MAX_PRIMES
  unsigned num_primes;

  // the actual primes, in decreasing order
  uint64_t p[ZZ_MAX_PRIMES];

  // corresponding values of pinv and pinvb
  uint64_t pinv[ZZ_MAX_PRIMES];
  uint64_t pinvb[ZZ_MAX_PRIMES];

  // -------- precomputed data for performing CRT

  // s[i] = -(2^64)^(i-1) / (p[0]*...*p[i-1]) mod p[i]
  // for 1 <= i < ZZ_MAX_PRIMES
  uint64_t s[ZZ_MAX_PRIMES];

  // spinv[i] = floor(2^64 * s[i] / p[i])
  // for 1 <= i < ZZ_MAX_PRIMES
  uint64_t spinv[ZZ_MAX_PRIMES];

  // u[i] = p[0] * ... * p[i-1], stored in i limbs,
  // for 1 <= i <= num_primes
  mp_limb_t u[ZZ_MAX_PRIMES + 1][ZZ_MAX_PRIMES];

  // uhalf[i] = floor(p[0] * ... * p[i-1] / 2), stored in i limbs,
  // for 1 <= i <= num_primes
  mp_limb_t uhalf[ZZ_MAX_PRIMES + 1][ZZ_MAX_PRIMES];

  // -------- precomputed data for performing FFTs

  fft62_mod_t fft62_mod[ZZ_MAX_PRIMES];

} zz_moduli_t;


// initialises moduli with num_primes primes
void zz_moduli_init(zz_moduli_t* moduli, unsigned num_primes);

// destroys moduli
void zz_moduli_clear(zz_moduli_t* moduli);


#endif
