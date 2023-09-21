/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_MPNFFT_MOD_H
#define ZZ_MPNFFT_MOD_H

#include "mpnfft.h"
#include "fermat.h"


// computes approximate inverse of D = {dp,dn}, with dp[dn-1] != 0,
// i.e. I = {ip,in} such that
//          0 <= B^(dn+in-1) - I*D < 3*D
// Additionally, if dn <= in + 1, the output satisfies the stronger inequality
//          0 <= B^(dn+in-1) - I*D < 2*D
// {ip,in} must not overlap {dp,dn}
void zz_invertappr(mp_limb_t* ip, size_t in, mp_limb_t* dp, size_t dn,
		   unsigned num_primes, zz_moduli_t* moduli, int threads);


typedef struct
{
  // number of FFT primes, must satisfy 1 <= num_primes <= ZZ_MAX_PRIMES
  unsigned num_primes;
  zz_moduli_t* moduli;

  // max_xn = maximum number of limbs in dividend
  // must have max_xn >= dn
  size_t max_xn;

  // {dp,dn} = divisor D
  // must have dn >= 1 and dp[dn-1] != 0
  mp_limb_t* dp;
  size_t dn;

  // transform of D
  zz_fermat_params_t params1;
  zz_fermat_poly_t DP;

  // {ip,in} = I = approximate inverse of D
  // satisfies 0 < B^(max_xn+1) - I*D < 4*D
  // in = max_xn - dn + 2
  // (note: we don't actually store I, just its transform)
  size_t in;

  // transform of I
  zz_mpnfft_params_t params2;
  zz_mpnfft_poly_t IP;
}
zz_mpnfft_mod_t;

// initialise "mod" for computing X mod D, where D = {dp,dn}, dn >= 1,
// dp[dn-1] != 0, and X has at most max_xn limbs
// must have max_xn >= dn
void zz_mpnfft_mod_init(zz_mpnfft_mod_t* mod, size_t max_xn,
			mp_limb_t* dp, size_t dn,
			unsigned num_primes, zz_moduli_t* moduli, int threads);

void zz_mpnfft_mod_clear(zz_mpnfft_mod_t* mod);


// computes X mod D, where X = sign * {xp,xn}, sign = 1 or -1.
// output written to {rp,dn}, exactly dn limbs.
// output in [0, D).
// must have 0 <= xn <= max_xn.
// ok for rp to alias xp.
void zz_mpnfft_mod_mod(zz_mpnfft_mod_t* mod, mp_limb_t* rp,
		       mp_limb_t* xp, size_t xn, int sign, int threads);


#endif
