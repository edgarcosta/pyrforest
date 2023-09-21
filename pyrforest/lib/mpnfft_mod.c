/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <assert.h>
#include "zzmisc.h"
#include "fermat.h"
#include "mpnfft_mod.h"
#include "zzmem.h"

// INVERTAPPR_THRESHOLD = threshold for switching from GMP division to
// our own newton iteration implementation
#if TEST
// for better test coverage
#define INVERTAPPR_THRESHOLD 4
#else
// completely made-up value
#define INVERTAPPR_THRESHOLD 1000
#endif
#if INVERTAPPR_THRESHOLD < 4
#error illegal value for INVERTAPPR_THRESHOLD
#endif


void zz_invertappr(mp_limb_t* ip, size_t in, mp_limb_t* dp, size_t dn,
		   unsigned num_primes, zz_moduli_t* moduli, int threads)
{
  assert(dn >= 1);
  assert(dp[dn - 1] != 0);
  assert(in >= 1);

  if (dn > in + 1)
    {
      // write D = D0 + D1*B^(dn-in-1) and divide B^(2*in) by D1 instead
      zz_invertappr(ip, in, dp + dn - in - 1, in + 1,
		    num_primes, moduli, threads);
      if (mpn_sub_1(ip, ip, in, 1) != 0)
	ip[0] = 0;

      // correctness proof:
      // the recursive call gives us J such that
      //           0 <= B^(2*in) - J*D1 < 2*D1
      // ==>       0 <= B^(dn+in-1) - J*D1*B^(dn-in-1) < 2*D1*B^(dn-in-1)
      // ==>   -J*D0 <= B^(dn+in-1) - J*D < 2*D
      // But J*D0 < B^in*B^(dn-in-1) = B^(dn-1) <= D, so
      //           0 <= B^(dn+in-1) - (J-1)*D < 3*D,
      // so we can take I = J - 1.

      // In the special case J == 0, we have
      //       B^(2*in) < 2*D1 <= 2*B^(in+1),
      // so must have in == 1, D1 > B^2/2, and hence D > B^dn/2,
      // which is equivalent to B^(dn+in-1) < 2*D, so we can simply take I = 0.
      return;
    }

  // now can assume dn <= in + 1

  // base case
  if (in <= INVERTAPPR_THRESHOLD)
    {
      // x = B^(dn+in-1) - 1
      size_t xn = dn + in - 1;
      mp_limb_t xp[xn];
      for (size_t j = 0; j < xn; j++)
	xp[j] = (mp_limb_signed_t) (-1);
      // todo: we should be using mpn_invertappr, not mpn_tdiv_qr.
      // Unfortunately the former is not public.
      mpn_tdiv_qr(ip, xp, 0, xp, xn, dp, dn);
      return;
    }

  // jn = ceil((in+3)/2)
  // jn < in since in >= 5
  size_t jn = in / 2 + 2;
  mp_limb_t* jp = ip + in - jn;

  // write D = G + H*B^(dn-hn), with 0 <= G < B^(dn-hn) and B^(hn-1) <= H < B^hn
  size_t hn = MIN(jn + 1, dn);
  mp_limb_t* hp = dp + dn - hn;

  // recursive call to get J = {jp,jn} such that
  //    0 <= B^(hn+jn-1) - J*H < 2*H
  // (bound is 2*H, not 3*H, because hn <= jn + 1)
  zz_invertappr(jp, jn, hp, hn, num_primes, moduli, threads);

  // now we have
  //         0 <= B^(dn+jn-1) - J*H*B^(dn-hn) < 2*H*B^(dn-hn)
  // ==>  -J*G <= B^(dn+jn-1) - J*D < 2*D

  // if hn == jn + 1 then J*G < B^(jn+dn-hn) <= B^(dn-1) <= D.
  // if hn == dn then G == 0, so still J*G < D.
  // In both cases we thus have
  //    0 <= B^(dn+jn-1) - (J-1)*D < 3*D.
  // Note that J cannot be zero; if it were, then we would have
  // B^(hn+jn-1) < 2*H < 2*B^hn, which is impossible since jn >= 2.
  // So we can replace J by J - 1, and then we will have
  //    0 <= B^(dn+jn-1) - J*D < 3*D.
  mpn_sub_1(jp, jp, jn, 1);

  // Let E = B^(dn+jn-1) - J*D, so 0 <= E < 3*D < 3*B^dn.

  // Write E = S + T*B^(dn-tn), with 0 <= S < B^(dn-tn) and 0 <= T < 3*B^tn.
  size_t tn = MIN(in - jn + 2, dn);    // tn >= 1

  // We will compute I = J*B^(in-jn) + floor(J*T/B^(2*jn-in+tn-1)).

  // Proof that this choice of I works: we have
  //      D*J*(B^(dn+jn-1) + E) = B^(2*dn+2*jn-2) - E^2

  // ==>  D*(J*B^(in-jn) + J*E/B^(2*jn-in+dn-1)) =
  //                              B^(dn+in-1) - E^2/B^(2*jn-in+dn-1)

  // ==>  D*(J*B^(in-jn) + J*T/B^(2*jn-in+tn-1)) =
  //                B^(dn+in-1) - E^2/B^(2*jn-in+dn-1) - D*J*S/B^(2*jn-in+dn-1)

  // ==>  D*I = B^(dn+in-1) - E^2/B^(2*jn-in+dn-1) - D*J*S/B^(2*jn-in+dn-1) - X
  // where 0 <= X < D.

  // to estimate the error terms we have
  //   E^2/B^(2*jn-in+dn-1) < 9*B^(2*dn)*B^(in-dn-2*jn+1)
  //                        <= 9*B^(dn+in-2*jn+1) <= 9*B^(dn-2) <= 9*D/B.
  // since 2*jn >= in + 3.

  // If tn == in - jn + 2 then
  //   J*S/B^(2*jn-in+dn-1) < B^jn*B^(dn-tn)*B^(in-dn-2*jn+1) <= B^(-1),
  // and if tn == dn then S = 0, so in both cases
  //   D*J*S/B^(2*jn-in+dn-1) < D/B.

  // Therefore the total error is at most
  //   0 <= B^(dn+in-1) - D*I < D*(1 + 9/B + 1/B) < 2*D.

  // Two possible strategies:
  //    reuse == 1: reuse fermat transform of J
  //    reuse == 0: recompute transform of J with different transform length
  // (the crossover is a first-order guess based on theoretical FFT cost)
  int reuse = (in < 3 * dn / 2);

  size_t yn = jn + tn + 1;

  // (1) compute E = B^(dn+jn-1) - J*D using fermat transform

  zz_fermat_params_t params1;
  zz_fermat_params_init(&params1, reuse ? MAX(dn + 1, yn) : (dn + 1),
			1, num_primes, moduli);

  size_t fermat_bufsize = zz_fermat_params_bufsize(&params1);
  mp_limb_t* buf = zz_malloc(MAX(fermat_bufsize, yn) * sizeof(mp_limb_t));
  mp_limb_t* yp = buf;    // later buffer will hold Y = J*T

  zz_fermat_transform_t transform;
  zz_fermat_transform_init(transform, &params1, buf);

  // todo: could do the following two transforms in parallel?
  zz_fermat_poly_t DP;
  zz_fermat_mpn_to_transform(transform, dp, dn, threads);
  zz_fermat_poly_init(DP, &params1);
  zz_fermat_transform_to_poly(DP, transform, -1, 64, 1, threads);   // -D
  zz_fermat_poly_fft(DP, DP, threads);

  zz_fermat_poly_t JP;
  zz_fermat_mpn_to_transform(transform, jp, jn, threads);
  zz_fermat_poly_init(JP, &params1);
  zz_fermat_transform_to_poly(JP, transform, 1, 0, 0, threads);
  zz_fermat_poly_fft(JP, JP, threads);

  zz_fermat_poly_mul(DP, DP, JP, 1, threads);
  zz_fermat_poly_ifft(DP, DP, threads);
  zz_fermat_poly_to_transform(transform, DP, threads);
  zz_fermat_poly_clear(DP);

  // now "transform" contains fermat transform of -J*D = E - B^(dn+jn-1).
  // add back B^(dn+jn-1) to each residue
  mp_limb_t* res = buf;
  for (unsigned i = 0; i < params1.w; i++)
    {
      // reduce B^(dn+jn-1) mod B^(M_i*t) + 1
      size_t len = params1.t << params1.e[i];    // M_i * t
      size_t k = (dn + jn - 1) % (2 * len);
      if (k < len)
	{
	  transform->hi[i] += mpn_add_1(res + k, res + k, len - k, 1);
	  if (transform->hi[i] == 2)     // overflow
	    transform->hi[i] = 1 - mpn_sub_1(res, res, len, 1);
	}
      else
	{
	  k -= len;
	  transform->hi[i] -= mpn_sub_1(res + k, res + k, len - k, 1);
	  if (transform->hi[i] == (mp_limb_signed_t) (-2))    // overflow
	    transform->hi[i] = mpn_add_1(res, res, len, 1) - 1;
	}
      res += len;
    }

  // put E into {buf, dn+1}
  zz_fermat_transform_to_mpn(buf, transform, threads);     // inplace
  // and hence T into {tp,tn+1}
  mp_limb_t* tp = buf + dn - tn;

  // (2) compute Y = J * T

  if (reuse)
    {
      // on this branch we use the same fermat transform parameters,
      // i.e. we can reuse the transform of J

      // move T to beginning of transform buffer
      mpn_copyi(buf, tp, tn + 1);
      tp = buf;

      zz_fermat_mpn_to_transform(transform, tp, tn + 1, threads);

      zz_fermat_poly_t TP;
      zz_fermat_poly_init(TP, &params1);
      zz_fermat_transform_to_poly(TP, transform, 1, 64, 1, threads);
      zz_fermat_poly_fft(TP, TP, threads);

      zz_fermat_poly_mul(TP, TP, JP, 1, threads);
      zz_fermat_poly_clear(JP);

      zz_fermat_poly_ifft(TP, TP, threads);
      zz_fermat_poly_to_transform(transform, TP, threads);
      zz_fermat_poly_clear(TP);

      // put Y into {yp,yn}
      zz_fermat_transform_to_mpn(yp, transform, threads);   // inplace
    }
  else
    {
      // on this branch we use a regular (non-fermat) transform,
      // so we need to transform J again

      zz_fermat_poly_clear(JP);

      zz_mpnfft_params_t params2;
      zz_mpnfft_params_init(&params2, 64 * yn, 1, num_primes, moduli);

      // todo: could do the following two transforms in parallel
      zz_mpnfft_poly_t TP;
      zz_mpnfft_poly_init(TP, &params2);
      zz_mpnfft_mpn_to_poly(TP, tp, tn + 1, 1, 64, 1, threads);
      zz_mpnfft_poly_fft(TP, TP, threads);

      zz_mpnfft_poly_t JP2;
      zz_mpnfft_poly_init(JP2, &params2);
      zz_mpnfft_mpn_to_poly(JP2, jp, jn, 1, 0, 0, threads);
      zz_mpnfft_poly_fft(JP2, JP2, threads);

      zz_mpnfft_poly_mul(JP2, JP2, TP, 1, threads);
      zz_mpnfft_poly_clear(TP);
      zz_mpnfft_poly_ifft(JP2, JP2, 0, threads);
      // put Y into {yp,yn}
      zz_mpnfft_poly_to_mpn(yp, yn, JP2, threads);
      zz_mpnfft_poly_clear(JP2);

      zz_mpnfft_params_clear(&params2);
    }

  // now we have Y in {yp,yn}

  // let Z = floor(Y / B^(2*jn-in+tn-1))
  // size of Z is yn - (2*jn - in + tn - 1) = in - jn + 2
  size_t zn = in - jn + 2;
  mp_limb_t* zp = yp + yn - zn;

  // (3) write J*B^(in-jn) + Z to output
  // todo: parallel copy
  mpn_copyi(ip, zp, zn - 2);
  if (mpn_add(jp, jp, jn, zp + zn - 2, 2))   // top two limbs of Z
    {
      // overflow; can just take I = B^in - 1
      for (size_t i = 0; i < in; i++)
      	ip[i] = (mp_limb_signed_t) (-1);
    }

  zz_fermat_transform_clear(transform);
  zz_fermat_params_clear(&params1);
  zz_free(buf, MAX(fermat_bufsize, yn) * sizeof(mp_limb_t));
}


void zz_mpnfft_mod_init(zz_mpnfft_mod_t* mod, size_t max_xn,
			mp_limb_t* dp, size_t dn,
			unsigned num_primes, zz_moduli_t* moduli, int threads)
{
  assert(dn >= 1);
  assert(dp[dn - 1] != 0);
  assert(max_xn >= dn);

  mod->max_xn = max_xn;
  mod->dn = dn;
  mod->num_primes = num_primes;
  mod->moduli = moduli;
  mod->in = max_xn - dn + 2;

  mod->dp = zz_malloc(dn * sizeof(mp_limb_t));
  mpn_copyi(mod->dp, dp, dn);

  // compute I satisfying 0 <= B^(max_xn+1) - I*D < 3*D
  mp_limb_t* ip = zz_malloc(mod->in * sizeof(mp_limb_t));
  zz_invertappr(ip, mod->in, mod->dp, dn, num_primes, moduli, threads);

  // I cannot be zero, since B^(max_xn+1) > 3*B^max_xn >= 3*D.
  // Replace I by I - 1, then we have 0 < B^(max_xn+1) - I*D < 4*D.
  mpn_sub_1(ip, ip, mod->in, 1);

  // compute transform of I
  zz_mpnfft_params_init(&mod->params2, 64 * 2 * (max_xn - dn + 2), 1,
			num_primes, moduli);
  zz_mpnfft_poly_init(mod->IP, &mod->params2);
  zz_mpnfft_mpn_to_poly(mod->IP, ip, mod->in, 1, 64, 1, threads);
  zz_free(ip, mod->in * sizeof(mp_limb_t));
  zz_mpnfft_poly_fft(mod->IP, mod->IP, threads);

  // compute transform of D
  zz_fermat_params_init(&mod->params1, dn + 1, 1, num_primes, moduli);
  zz_fermat_transform_t DT;
  zz_fermat_transform_init(DT, &mod->params1, NULL);
  zz_fermat_mpn_to_transform(DT, mod->dp, dn, threads);
  zz_fermat_poly_init(mod->DP, &mod->params1);
  zz_fermat_transform_to_poly(mod->DP, DT, 1, 64, 1, threads);
  zz_fermat_transform_clear(DT);
  zz_fermat_poly_fft(mod->DP, mod->DP, threads);
}


void zz_mpnfft_mod_clear(zz_mpnfft_mod_t* mod)
{
  zz_fermat_poly_clear(mod->DP);
  zz_fermat_params_clear(&mod->params1);
  zz_mpnfft_poly_clear(mod->IP);
  zz_mpnfft_params_clear(&mod->params2);
  zz_free(mod->dp, mod->dn * sizeof(mp_limb_t));
}


void zz_mpnfft_mod_mod(zz_mpnfft_mod_t* mod, mp_limb_t* rp,
		       mp_limb_t* xp, size_t xn, int sign, int threads)
{
  assert(xn <= mod->max_xn);

  size_t max_xn = mod->max_xn;
  size_t dn = mod->dn;

  if (xn < dn || dn <= 2)
    {
      if (xn < dn)
	{
	  // easy case: dividend is smaller than divisor.
	  if (rp != xp && xn != 0)
	    mpn_copyi(rp, xp, xn);
	  mpn_zero(rp + xn, dn - xn);
	}
      else
	{
	  // special case: very small divisor. Let GMP do it.
	  mp_limb_t* qp = zz_malloc((xn - dn + 1) * sizeof(mp_limb_t));
	  mpn_tdiv_qr(qp, rp, 0, xp, xn, mod->dp, dn);
	  zz_free(qp, (xn - dn + 1) * sizeof(mp_limb_t));
	}

      // replace result by D - {rp,dn} if necessary
      if (sign == -1)
	{
	  int zero = 1;
	  for (size_t i = 0; i < dn; i++)
	    if (rp[i] != 0)
	      {
		zero = 0;
		break;
	      }

	  if (!zero)
	    mpn_sub_n(rp, mod->dp, rp, dn);
	}

      return;
    }

  // H = {hp,hn} = floor(|X| / B^(dn-2))
  mp_limb_t* hp = xp + dn - 2;
  size_t hn = xn - dn + 2;

  size_t bufsize = zz_fermat_params_bufsize(&mod->params1);
  size_t yn = hn + mod->in;
  size_t qn = yn - max_xn + dn - 3;
  mp_limb_t* buf = zz_malloc((yn - qn + MAX(bufsize, qn + 1)) * sizeof(mp_limb_t));
  mp_limb_t* yp = buf;
  mp_limb_t* qp = buf + yn - qn;

  // compute Y = H * I and then Q = floor(H * I / B^(max_xn-dn+3))
  zz_mpnfft_poly_t HP;
  zz_mpnfft_poly_init(HP, &mod->params2);
  zz_mpnfft_mpn_to_poly(HP, hp, hn, 1, 0, 0, threads);
  zz_mpnfft_poly_fft(HP, HP, threads);
  zz_mpnfft_poly_mul(HP, HP, mod->IP, 1, threads);
  zz_mpnfft_poly_ifft(HP, HP, 0, threads);
  zz_mpnfft_poly_to_mpn(yp, yn, HP, threads);
  // now Q is in {qp,qn}
  zz_mpnfft_poly_clear(HP);

  // Claim: at this stage 0 < |X| - Q*D < 2D.
  // Proof: we have
  //    0 < B^(max_xn+1) - D*I < 4*D,
  //    0 <= |X| - H*B^(dn-2) < B^(dn-2) <= D/B,
  //    0 <= H*I/B^(max_xn-dn+3) - Q < 1.
  // Multiplying these by respectively H/B^(max_xn-dn+3) < 1/B, 1, D,
  // and adding together, yields
  //    0 < |X| - Q*D < 4D/B + D/B + D < 2D.

  // Now compute
  //    |X| - Q*D  if sign == 1, or
  //    (Q+1)*D - |X|  if sign == -1,
  // modulo the fermat modulus.
  // In the first case the result is in (0, 2*D) (and most likely in (0, D))
  // In the second case the result is in (-D, D) (and most likely in (0, D))

  zz_fermat_transform_t QT;
  if (sign == -1)
    {
      // replace Q by Q + 1
      // (possibly lengthen Q buffer if overflow occurs)
      if (mpn_add_1(qp, qp, qn, 1))
	qp[qn++] = 1;
    }

  zz_fermat_transform_init(QT, &mod->params1, qp);
  zz_fermat_mpn_to_transform(QT, qp, qn, threads);   // inplace
  zz_fermat_poly_t QP;
  zz_fermat_poly_init(QP, &mod->params1);
  zz_fermat_transform_to_poly(QP, QT, 1, 0, 0, threads);
  zz_fermat_poly_fft(QP, QP, threads);
  zz_fermat_poly_mul(QP, QP, mod->DP, 1, threads);
  zz_fermat_poly_ifft(QP, QP, threads);
  zz_fermat_poly_to_transform(QT, QP, threads);
  zz_fermat_poly_clear(QP);

  zz_fermat_transform_t XT;
  zz_fermat_transform_init(XT, &mod->params1, NULL);
  zz_fermat_mpn_to_transform(XT, xp, xn, threads);
  if (sign == 1)
    zz_fermat_transform_sub(XT, XT, QT, threads);
  else
    zz_fermat_transform_sub(XT, QT, XT, threads);
  zz_fermat_transform_clear(QT);
  zz_free(buf, (yn - qn + MAX(bufsize, qn + 1)) * sizeof(mp_limb_t));
  zz_fermat_transform_to_mpn(XT->data, XT, threads);

  // move result to rp, with final unlikely reduction if necessary
  if (sign == 1)
    {
      // if result >= D, subtract D
      if (XT->data[dn] != 0 || mpn_cmp(XT->data, mod->dp, mod->dn) >= 0)
	// must subtract D once
	mpn_sub_n(rp, XT->data, mod->dp, dn);
      else
	// result is correct already
	mpn_copyi(rp, XT->data, dn);
    }
  else
    {
      // if result < 0, add D (modulo fermat modulus)
      size_t Fsize = mod->params1.M * mod->params1.t;
      if (XT->data[Fsize] || XT->data[Fsize - 1] >> 63)
	{
	  zz_fermat_add_F(XT->data, -1, &mod->params1);
	  mpn_add_n(rp, XT->data, mod->dp, dn);
	}
      else
	// result is correct already
	mpn_copyi(rp, XT->data, dn);
    }

  zz_fermat_transform_clear(XT);
}
