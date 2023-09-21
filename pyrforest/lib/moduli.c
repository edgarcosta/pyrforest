/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <assert.h>
#include "fft62/mod62.h"
#include "mpzfft_moduli.h"


uint64_t global_p[ZZ_MAX_PRIMES] =
  {FFT62_PRIME1, FFT62_PRIME2, FFT62_PRIME3, FFT62_PRIME4,
   FFT62_PRIME5, FFT62_PRIME6, FFT62_PRIME7, FFT62_PRIME8};


void zz_moduli_init(zz_moduli_t* moduli, unsigned num_primes)
{
  assert(1 <= num_primes && num_primes <= ZZ_MAX_PRIMES);
  moduli->num_primes = num_primes;

  for (unsigned i = 0; i < num_primes; i++)
    assert(mod62_valid(global_p[i]));
  for (unsigned i = 1; i < num_primes; i++)
    assert(global_p[i] < global_p[i-1]);

  for (unsigned i = 0; i < num_primes; i++)
    {
      uint64_t p = global_p[i];
      moduli->p[i] = p;
      moduli->pinvb[i] = mod62_pinvb(p);
      moduli->pinv[i] = mod62_pinv(p);
      fft62_mod_init(&moduli->fft62_mod[i], p);
    }

  for (unsigned i = 1; i <= num_primes; i++)
    {
      if (i == 1)
	moduli->u[1][0] = global_p[0];
      else
	moduli->u[i][i-1] =
	  mpn_mul_1(moduli->u[i], moduli->u[i-1], i-1, global_p[i-1]);

      mpn_rshift(moduli->uhalf[i], moduli->u[i], i, 1);
    }

  for (unsigned i = 1; i < num_primes; i++)
    {
      uint64_t p = moduli->p[i];
      uint64_t pinv = moduli->pinv[i];
      uint64_t s = mpn_mod_1(moduli->u[i], i, p);
      s = p - mod62_pow_pinv(s, p - 2, p, pinv);
      for (unsigned j = 0; j < i - 1; j++)
	s = mod62_mul_pinv(s, (uint64_t) (-p), p, pinv);
      moduli->s[i] = s;
      moduli->spinv[i] = mod62_ypinv(s, p, pinv);
    }
}


void zz_moduli_clear(zz_moduli_t* moduli)
{
  for (unsigned i = 0; i < moduli->num_primes; i++)
    fft62_mod_clear(&moduli->fft62_mod[i]);
}
