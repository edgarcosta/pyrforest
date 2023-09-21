/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <stddef.h>
#include "fft62/mod62.h"
#include "mpzfft_moduli.h"


void zz_crt(mp_limb_t* rp, uint64_t** up, size_t un,
	    zz_moduli_t* moduli, unsigned num_primes)
{
  if (num_primes == 1)
    {
      // only one prime, simply reduce data mod p
      uint64_t p = moduli->p[0];
      uint64_t* src = up[0];
      for (ptrdiff_t j = un - 1; j >= 0; j--)
	{
	  uint64_t x = mod62_reduce4(src[j], p);
	  rp[j] = x - ((x > p/2) ? p : 0);
	}
      return;
    }

  // first pass: CRT data from first two primes
  {
    uint64_t* up0 = up[0];
    uint64_t* up1 = up[1];
    mp_limb_t* dest = rp + (un - 1) * num_primes;
    uint64_t p0 = moduli->p[0];
    uint64_t p1 = moduli->p[1];
    uint64_t s = moduli->s[1];
    uint64_t spinv = moduli->spinv[1];

    for (ptrdiff_t j = un - 1; j >= 0; j--, dest -= num_primes)
      {
	uint64_t x = mod62_reduce4(up0[j], p0);       // [0, p0)
	uint64_t y = up1[j];                          // [0, 4*p1)
	uint64_t r = x - y + ((x < y) ? (4*p1) : 0);  // [0, B)
	r = mod62_mul_ypinv(r, s, spinv, p1);         // r = (y - x)/p0 mod p1
	uint64_t s1, s0;
	MUL_WIDE(s1, s0, r, p0);
	s0 += x;
	dest[0] = s0;
	dest[1] = s1 + (s0 < x);    // carry
      }
  }

  // remaining primes
  for (unsigned i = 2; i < num_primes; i++)
    {
      uint64_t* src = up[i];
      mp_limb_t* dest = rp + (un - 1) * num_primes;
      uint64_t p = moduli->p[i];
      uint64_t pinvb = moduli->pinvb[i];
      mp_limb_t* u = moduli->u[i];
      uint64_t s = moduli->s[i];
      uint64_t spinv = moduli->spinv[i];

      for (ptrdiff_t j = un - 1; j >= 0; j--, dest -= num_primes)
	{
	  // Let X be the previous residue, stored in the first i limbs of dest,
	  // i.e. X = up[k][j] mod p[k] for 0 <= k < i, and
	  // 0 <= X < p[0]*...*p[i-1].

	  // Compute r = (X - up[i][j]) / B^(i-1) mod p[i], in [0, 4*p[i]).
	  uint64_t x = dest[0];     // [0, B)
	  uint64_t y = src[j];      // [0, 4p)
	  uint64_t r = x - y + ((x < y) ? (4*p) : 0);
	  for (unsigned k = 1; k < i; k++)
	    {
	      uint64_t h = MUL_HI(r * pinvb, p);      // [0, p)
	      uint64_t x = dest[k];                   // [0, B)
	      r = x - h + ((x < h) ? p : 0);
	    }

	  // r := (up[i][j] - X) / (p[0]*...*p[i-1]) mod p[i], in [0, p[i])
	  r = mod62_mul_ypinv(r, s, spinv, p);

	  // add r*p[0]*...*p[i-1] to X
	  dest[i] = mpn_addmul_1(dest, u, i, r);
	}
    }

  // convert to negative residues if required
  mp_limb_t* dest = rp + (un - 1) * num_primes;
  mp_limb_t* u = moduli->u[num_primes];
  mp_limb_t* uhalf = moduli->uhalf[num_primes];
  for (ptrdiff_t j = un - 1; j >= 0; j--, dest -= num_primes)
    if (mpn_cmp(dest, uhalf, num_primes) > 0)
      mpn_sub_n(dest, dest, u, num_primes);
}
