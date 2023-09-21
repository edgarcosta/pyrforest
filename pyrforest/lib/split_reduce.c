/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <stddef.h>
#include "fft62/mod62.h"
#include "mpzfft_moduli.h"
#include "split.h"
#include "reduce.h"
#include "zzmisc.h"
#include "zzmem.h"


void zz_split_reduce(uint64_t** rp, size_t rn, int sign, int lgS,
		     zz_moduli_t* moduli, unsigned num_primes,
		     const mp_limb_t* up, size_t un, unsigned b,
		     int threads)
{
  // limbs per split coefficient
  unsigned t = (b + 63) / 64;

  // chunk size, want everything to safely fit into L1
  // (must be divisible by 64)
  const size_t size = 512;

  // c[j] = sign * 2^lgS * B^(t-1) mod p[j]
  uint64_t c[ZZ_MAX_PRIMES];
  uint64_t cpinv[ZZ_MAX_PRIMES];
  for (unsigned j = 0; j < num_primes; j++)
    {
      uint64_t p = moduli->p[j];
      uint64_t pinv = moduli->pinv[j];
      uint64_t x = mod62_2exp(lgS + 64*(t - 1), p, pinv);
      x = (sign == 1) ? x : (p - x);
      c[j] = x;
      cpinv[j] = mod62_ypinv(x, p, pinv);
    }

//#pragma omp parallel num_threads(threads)
  {
    mp_limb_t* temp = zz_malloc(t * size * sizeof(mp_limb_t));

  // process {rp,rn} in chunks
//#pragma omp for schedule(static)
    for (ptrdiff_t i0 = 0; i0 < rn; i0 += size)
      {
	// this iteration handles outputs i0 <= i < i1
	size_t i1 = MIN(i0 + size, rn);

	// split into integer coefficients
	ptrdiff_t offset = (i0 / 64) * b;
	if (un > offset)
	  zz_split(temp, i1 - i0, up + offset, un - offset, b);
	else
	  mpn_zero(temp, (i1 - i0) * t);

	// reduce modulo each prime
	for (unsigned j = 0; j < num_primes; j++)
	  zz_reduce(rp[j] + i0, temp, i1 - i0, t, c[j], cpinv[j],
		    moduli->p[j], moduli->pinvb[j]);
      }

    zz_free(temp, t * size * sizeof(mp_limb_t));
  }
}
