/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_REDUCE_H
#define ZZ_REDUCE_H

#include <stdint.h>
#include <gmp.h>


/*
  INPUT: t >= 1
         {up, t*un}, un unsigned integers of t limbs, say d[0], ... d[un-1]
         c in [0, p)
         cpinv = floor(y * B / p)
         pinvb = 1 / p mod B
  OUTPUT: rp[i] = c * d[i] / B^(t-1) mod p, in [0, 2p), for 0 <= i < un.
*/
void zz_reduce(uint64_t* rp, const mp_limb_t* up, size_t un, unsigned t,
	       uint64_t c, uint64_t cpinv, uint64_t p, uint64_t pinvb);


#endif
