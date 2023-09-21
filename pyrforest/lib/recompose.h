/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_RECOMPOSE_H
#define ZZ_RECOMPOSE_H

#include <gmp.h>


/*
  INPUT: rn >= 1
         un >= 0
	 r0 >= 0
         r >= 1
         s >= r / 64
  OUTPUT:
  Interprets input as a sequence of un signed integers, s limbs each, stored
  consecutively in up, in 2's complement, say d[0], ... d[un-1].
  Output is 2^r0 * (d[0] + d[1]*2^r + ... + d[un-1]*2^((un-1)r)) mod B^rn,
  stored at rp.
*/
void zz_recompose(mp_limb_t* rp, size_t rn, mp_limb_t* up, size_t un,
		  unsigned s, unsigned r, unsigned r0);

#endif
