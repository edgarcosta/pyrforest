/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <stddef.h>
#include "zzmisc.h"


void zz_split(mp_limb_t* rp, size_t rn, const mp_limb_t* up, size_t un,
	      unsigned b)
{
  if (b % 64 == 0)
    {
      // split along limb boundaries
      size_t k = (b / 64) * rn;
      size_t d = MIN(k, un);
      mpn_copyi(rp, up, d);
      if (k > d)
	mpn_zero(rp + d, k - d);
    }
  else
    {
      // t = number of limbs to write for each output coefficient
      size_t t = b / 64 + 1;
      b %= 64;

      // for last limb of each output coefficient
      mp_limb_t mask = ((mp_limb_t) 1 << b) - 1;

      // each iteration reads from bit index "bit" at "up"
      unsigned bit = 0;

      // main splitting loop
      while (rn > 0 && (ptrdiff_t) un > t)
	{
	  if (bit != 0)
	    {
	      mpn_rshift(rp, up, t, bit);
	      rp[t - 1] |= up[t] << (64 - bit);
	    }
	  else
	    mpn_copyi(rp, up, t);
	  rp[t - 1] &= mask;

	  rp += t;
	  rn--;

	  bit += b;
	  size_t m = (t - 1) + (bit / 64);
	  up += m;
	  un -= m;
	  bit %= 64;
	}

      // last portion of input
      while (rn > 0 && (ptrdiff_t) un > 0)
	{
	  if (bit != 0)
	    mpn_rshift(rp, up, un, bit);
	  else
	    mpn_copyi(rp, up, un);
	  if (un < t)
	    mpn_zero(rp + un, t - un);
	  rp[t - 1] &= mask;

	  rp += t;
	  rn--;

	  bit += b;
	  size_t m = (t - 1) + (bit / 64);
	  up += m;
	  un -= m;
	  bit %= 64;
	}

      // last portion of output
      if (rn > 0)
	mpn_zero(rp, rn * t);
    }
}
