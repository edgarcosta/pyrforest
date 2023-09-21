/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_CRT_H
#define ZZ_CRT_H

#include "mpzfft_moduli.h"


/*
  INPUT: un >= 1
         num_primes >= 1, number of CRT primes
         up = arrays up[0], up[1], ..., up[num_primes-1], each of length un
         moduli = moduli object with at least num_primes primes
  OUTPUT:
  Entries of up[i] are interpreted mod p[i] (= moduli->p[i]), in [0, 4p).
  For each 0 <= j < un, computes CRT of up[0][j], ..., up[num_primes-1][j],
  writes result (num_primes limbs) starting at rp[j*num_primes], i.e. total
  output is num_primes * un limbs.
  Each output coefficient is reduced into the canonical interval -P/2 < x < P/2,
  where P = p[0] * ... * p[num_primes-1], and written in 2's complement.
*/
void zz_crt(mp_limb_t* rp, uint64_t** up, size_t un,
	    zz_moduli_t* moduli, unsigned num_primes);


#endif
