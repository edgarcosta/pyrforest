/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_SPLIT_H
#define ZZ_SPLIT_H


#include <gmp.h>


/*
  INPUT: un >= 1
         b >= 1
         rn >= 1
  OUTPUT:
  Let
     {up,un} = d[0] + 2^b*d[1] + ... + 2^(b*(rn-1))*d[rn-1]  mod 2^(b*rn),
  where each 0 <= d[i] < 2^b. Output is {rp,t*rn} = d[0], ..., d[rn-1], where
  t = ceil(b / 64), i.e. t limbs per output coefficient, each unsigned.
*/
void zz_split(mp_limb_t* rp, size_t rn, const mp_limb_t* up, size_t un,
	      unsigned b);

#endif
