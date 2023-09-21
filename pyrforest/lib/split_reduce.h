/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_SPLIT_REDUCE_H
#define ZZ_SPLIT_REDUCE_H

#include <gmp.h>
#include "mpzfft_moduli.h"


/*
  Combined operation of zz_split() and zz_reduce() for num_primes primes.

  More precisely, input and coefficient splitting is described by the parameters
  up, un, b, as in zz_split(). Let d[i] be the i-th split coefficient.

  For each 0 <= j < num_primes, writes sign * d[i] * 2^lgS mod p[j] to rp[j][i]
  for 0 <= i < rn, where sign = -1 or 1.

  Outputs in [0, 2p).

  Data is processed in blocks for better locality.

  threads = maximum number of threads to use.
*/
void zz_split_reduce(uint64_t** rp, size_t rn, int sign, int lgS,
		     zz_moduli_t* moduli, unsigned num_primes,
		     const mp_limb_t* up, size_t un, unsigned b,
		     int threads);


#endif
