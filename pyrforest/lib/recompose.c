/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <assert.h>
#include <stddef.h>
#include "zzmisc.h"


void zz_recompose(mp_limb_t* rp, size_t rn, mp_limb_t* up, size_t un,
		  unsigned s, unsigned r, unsigned r0)
{
  // reduce to case r0 < 64
  while (r0 >= 64)
    {
      r0 -= 64;
      *rp++ = 0;
      if (--rn == 0)
	return;
    }

  mp_limb_t temp[s + 1];

  // on each iteration of main loop, dst increases by either u or u + 1
  size_t u = r / 64;
  assert(s >= u);

  ptrdiff_t src = 0;
  ptrdiff_t dst = 0;
  unsigned dst_bit = r0;

  // at the beginning of each iteration, the first dst + s - u + 1 limbs of
  // rp contain valid data
  mpn_zero(rp, MIN(s - u + 1, rn));
  mp_limb_t rp_cy = 0;     // implicit borrow at index dst + s - u + 1
  while (src < un && dst < (ptrdiff_t) (rn - s - 1))
    {
      mp_limb_t up_cy = up[s-1] >> 63;   // borrow for input coefficient
      if (dst_bit != 0)
	temp[s] = mpn_lshift(temp, up, s, dst_bit) + ((-up_cy) << dst_bit);
      else
	{
	  mpn_copyi(temp, up, s);
	  temp[s] = -up_cy;
	}

      mp_limb_t cy = mpn_add_n(rp + dst, rp + dst, temp, s - u + 1);
      if (u > 0)
	{
	  cy = mpn_add_1(rp + dst + s - u + 1, temp + s - u + 1, u, cy);
	  cy -= mpn_sub_1(rp + dst + s - u + 1, rp + dst + s - u + 1, u, rp_cy);
	}
      else
	cy -= rp_cy;
      cy -= up_cy;
      rp[dst + s + 1] = cy;
      rp_cy = -cy;
      assert(rp_cy == 0 || rp_cy == 1);

      dst_bit += r;
      dst += dst_bit / 64;
      dst_bit &= 63;
      src++;
      up += s;
    }

  // same loop as above, but don't write beyond end of {rp,rn}
  while (src < un && dst < rn)
    {
      mp_limb_t up_cy = up[s-1] >> 63;   // borrow for input coefficient
      if (dst_bit != 0)
	temp[s] = mpn_lshift(temp, up, s, dst_bit) + ((-up_cy) << dst_bit);
      else
	{
	  mpn_copyi(temp, up, s);
	  temp[s] = -up_cy;
	}

      mp_limb_t cy = mpn_add_n(rp + dst, rp + dst, temp,
			       MIN(s - u + 1, rn - dst));
      if (u > 0)
	{
	  if (dst + s - u + 1 < rn)
	    {
	      size_t k = MIN(u, rn - (dst + s - u + 1));
	      cy = mpn_add_1(rp + dst + s - u + 1, temp + s - u + 1, k, cy);
	      cy -= mpn_sub_1(rp + dst + s - u + 1, rp + dst + s - u + 1,
			      k, rp_cy);
	    }
	}
      else
	cy -= rp_cy;
      cy -= up_cy;
      rp_cy = -cy;

      dst_bit += r;
      dst += dst_bit / 64;
      dst_bit &= 63;
      src++;
      up += s;
    }

  // propagate final borrow to rest of output
  for (ptrdiff_t i = dst + s - u + 1; i < (ptrdiff_t) rn; i++)
    rp[i] = -rp_cy;
}
