/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_CRT_RECOMPOSE_H
#define ZZ_CRT_RECOMPOSE_H

#include "mpzfft_moduli.h"


/*
  Combined operation of zz_crt() and zz_recompose() for num_primes primes.

  Input is described by parameters up, un, moduli, num_primes as in zz_crt().

  Output is described by parameters rp, rn, r as in zz_recompose()
  (with r0 == 0).

  Must have rn >= 1, num_primes >= r / 64.

  Data is processed in blocks for better locality.

  threads = maximum number of threads to use.
*/
void zz_crt_recompose(mp_limb_t* rp, size_t rn, unsigned r,
		      uint64_t** up, size_t un,
		      zz_moduli_t* moduli, unsigned num_primes, int threads);

#endif
