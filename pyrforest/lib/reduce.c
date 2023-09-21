/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include "fft62/mod62.h"
#include <gmp.h>


void zz_reduce(uint64_t* rp, const mp_limb_t* up, size_t un, unsigned t,
	       uint64_t c, uint64_t cpinv, uint64_t p, uint64_t pinvb)
{
  while (un > 0)
    {
      uint64_t acc = up[0];

      for (size_t i = 1; i < t; i++)
	{
	  // acc = (up[0] + ... + up[i]*B^i) / B^i mod p,  0 <= acc < B
	  acc = MUL_HI(acc * pinvb, p);
	  acc = (up[i] - acc) + ((acc > up[i]) ? p : 0);
	}

      *rp++ = mod62_mul_ypinv_lazy(acc, c, cpinv, p);
      up += t;
      un--;
    }
}
