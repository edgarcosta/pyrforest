/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <assert.h>
#include "mpzfft.h"


void mpzfft_fft(mpzfft_t rop, mpz_t op, int threads)
{
  int sign = mpz_sgn(op);

  if (sign == 0)
    {
      rop->size = 0;
    }
  else
    {
      zz_mpnfft_mpn_to_poly(rop, op->_mp_d, mpz_size(op), sign, 0, 0, threads);
      zz_mpnfft_poly_fft(rop, rop, threads);
    }
}


void mpzfft_ifft(mpz_t rop, mpzfft_t op, int threads)
{
  if (op->size == 0)
    {
      mpz_set_ui(rop, 0);
      return;
    }

  zz_mpnfft_poly_ifft(op, op, 1, threads);

  // estimate of output size
  size_t n =
    ((op->size - 1) * op->params->r + 62 * op->params->num_primes + 2) / 64 + 1;
  mpz_realloc(rop, n);

  zz_mpnfft_poly_to_mpn(rop->_mp_d, n, op, threads);

  // complement if high bit set
  int neg = 0;
  if (rop->_mp_d[n - 1] >> 63)
    {
      // todo: parallelise
      mpn_neg(rop->_mp_d, rop->_mp_d, n);
      neg = 1;
    }
  // normalise
  while (n > 0 && rop->_mp_d[n - 1] == 0)
    n--;
  rop->_mp_size = neg ? -n : n;
}



void mpzfft_mod_init(mpzfft_mod_t* mod, size_t bits, mpz_t d,
		     unsigned num_primes, zz_moduli_t* moduli, int threads)
{
  assert(mpz_sgn(d) == 1);
  assert(bits >= 1);

  zz_mpnfft_mod_init(mod, (bits + 63) / 64, d->_mp_d, mpz_size(d),
		     num_primes, moduli, threads);
}



void mpzfft_mod_mod(mpzfft_mod_t* mod, mpz_t rop, mpz_t op, int threads)
{
  size_t dn = mod->dn;
  int sign = mpz_sgn(op);

  if (sign == 0)
    {
      mpz_set_ui(rop, 0);
      return;
    }

  if (rop->_mp_alloc < dn)
    {
      // not enough room in rop; need to grow it
      if (rop != op)
	{
	  // not inplace, can just destroy current value
	  mpz_set_ui(rop, 0);
	}
      mpz_realloc(rop, dn);
    }

  zz_mpnfft_mod_mod(mod, rop->_mp_d, op->_mp_d, mpz_size(op), sign, threads);

  // normalise
  size_t n = mod->dn;
  while (n > 0 && rop->_mp_d[n - 1] == 0)
    n--;
  rop->_mp_size = n;
}
