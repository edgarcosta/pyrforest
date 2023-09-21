/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <assert.h>
#include <math.h>
#include <string.h>
#include "zzmisc.h"
#include "fermat.h"
#include "split_reduce.h"
#include "crt_recompose.h"
#include "fft62/fft62.h"
#include "fft62/mod62.h"
#include "zzmem.h"


void zz_fermat_params_init(zz_fermat_params_t* params, size_t n, size_t terms,
			   unsigned num_primes, zz_moduli_t* moduli)
{
  params->moduli = moduli;
  params->num_primes = num_primes;
  params->terms = terms;

  // bound on polynomial coefficients
  // todo: this will overflow IEEE double after about 16 primes or so
  double bound = 0.4999;
  for (unsigned i = 0; i < num_primes; i++)
    bound *= (double) moduli->p[i];

  // try decreasing values of r until everything fits
  // todo: better initial estimate of r
  for (unsigned r = 62 * num_primes / 2; r > 0; r--)
    {
      // len = total number of polynomial coefficients, len >= 64 * n / r
      size_t len = ceil(64.0 * n / r);

      // choose M and K = 2^lgK with M not too large, M * K >= len,
      // r*K divisible by 64
      unsigned lgK = 0;
      while (((size_t) r << lgK) % 64)
	lgK++;
      // (threshold in next line controls granularity vs overhead of fermat
      // transform)
      while ((len >> lgK) > 64)
	lgK++;

      size_t K = (size_t) 1 << lgK;
      size_t t = ((size_t) r << lgK) / 64;
      size_t M = (len + K - 1) >> lgK;     // ceil(len / K)
      // b = number of bits in M
      unsigned b = 0;
      for (size_t temp = M; temp > 0; b++, temp >>= 1);
      assert(b <= ZZ_FERMAT_MAX_WEIGHT);

      // check whether 2^(2r) * terms * K * M_0 fits
      if (ldexp((double) terms, 2*r + lgK + b - 1) < bound)
	{
	  // yes, we're done
	  params->r = r;
	  params->lgK = lgK;
	  params->K = K;
	  params->t = t;
	  params->M = M;

	  // compute w and e[]
	  unsigned w = 0;
	  for (int i = b - 1; i >= 0; i--)
	    {
	      if ((M >> i) & 1)
		params->e[w++] = i;
	    }
	  params->w = w;

	  return;
	}
    }

  abort();
}



void zz_fermat_params_F_mpn(zz_fermat_params_t* params, mp_limb_t* rp)
{
  size_t t = params->t;
  size_t M = params->M;

  mpn_zero(rp, t * M);

  for (size_t i = 0; i <= M; i++)
    rp[i * t] = ((i & M) == i);
}


void zz_fermat_params_F_mpz(zz_fermat_params_t* params, mpz_t rop)
{
  size_t t = params->t;
  size_t M = params->M;

  mpz_realloc(rop, t * M + 1);
  zz_fermat_params_F_mpn(params, rop->_mp_d);
  rop->_mp_size = t * M + 1;
}



void zz_fermat_transform_init(zz_fermat_transform_t T,
			      zz_fermat_params_t* params, mp_limb_t* data)
{
  T->params = params;

  if (data == NULL)
    {
      size_t size = zz_fermat_params_bufsize(params);
      T->data = (mp_limb_t*) zz_malloc(size * sizeof(mp_limb_t));
      T->own = 1;
    }
  else
    {
      T->data = data;
      T->own = 0;
    }
}


void zz_fermat_transform_clear(zz_fermat_transform_t T)
{
  if (T->own)
    zz_free(T->data, zz_fermat_params_bufsize(T->params) * sizeof(mp_limb_t));
}


/*
  fwd_impl_* are building blocks for performing forward fermat transform.

  Input is a polynomial F(x) in Z[x]. Coefficients are chunks of size "chunk"
  limbs, as portions of {up,un}, starting at intervals of "skip" limbs. Here
  un == 0 is allowed. Any part of a coefficient not lying in {up,un} is assumed
  zero.

  Output is the sequence of residues
         G_0(x) = F(x) mod x^{M_0} + 1,
               ...
     G_{w-1}(x) = F(x) mod x^{M_{w-1}} + 1.

  The M output coefficients are stored as a main part of size "chunk" limbs,
  plus a signed carry limb. The main parts are written to
     {rp, chunk}
     {rp + skip, chunk}
     ...
     {rp + (M-1)*skip, chunk}.
  The carry limbs are written to
     rcy[0], rcy[skip_cy], ..., rcy[(M-1)*skip_cy].

  Inplace operation (i.e. up == rp) is allowed.

  scratch = scratch space of size ceil(log2(M)) * chunk.

  fwd_impl_2() additionally computes F(1), writes it to {sp, chunk}, and
  carry limb to scy. It also allows M == 0, i.e. compute only F(1).
*/

static void
fwd_impl_1(mp_limb_t* rp, mp_limb_signed_t* rcy,
	   mp_limb_t* up, size_t un,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch);

static void
fwd_impl_2(mp_limb_t* rp, mp_limb_signed_t* rcy,
	   mp_limb_t* sp, mp_limb_signed_t* scy,
	   mp_limb_t* up, size_t un,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch)
{
  if (un == 0)
    {
      // empty input, so all output coefficients are zero
      for (size_t i = 0; i < M; i++)
      	mpn_zero(rp + i * skip, chunk);
      for (size_t i = 0; i < M; i++)
	rcy[i * skip_cy] = 0;

      mpn_zero(sp, chunk);
      *scy = 0;

      return;
    }

  // wp = where to put F(-1) (if M odd)
  mp_limb_t* wp = rp + (M-1) * skip;
  mp_limb_signed_t* wcy = rcy + (M-1) * skip_cy;

  if (un <= skip)
    {
      // exactly one input coefficient

      if (M % 2)
	{
	  // F(-1) = u_0
	  if (wp != up)
	    mpn_copyi(wp, up, MIN(un, chunk));
	  if (chunk > un)
	    mpn_zero(wp + un, chunk - un);
	  *wcy = 0;
	}

      // F(1) = u_0
      if (sp != up)
	mpn_copyi(sp, up, MIN(un, chunk));
      if (chunk > un)
	mpn_zero(sp + un, chunk - un);
      *scy = 0;

      // recurse into even part
      fwd_impl_1(rp, rcy, up, un, M / 2, 2 * skip, 2 * skip_cy,
		 chunk, scratch);

      // recurse into odd part
      fwd_impl_1(rp + skip, rcy + skip_cy, up + skip, MAX(un, skip) - skip,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
  else if (un <= 2*skip)
    {
      // exactly two input coefficients

      if (M == 1)
	{
	  // F(-1) = u_0 - u_1, in temporary
	  mp_limb_t* tp = scratch;
	  mp_limb_signed_t tcy = -mpn_sub(tp, up, chunk,
					  up + skip, MIN(un - skip, chunk));

	  // F(1) = u_0 + u_1
	  *scy = mpn_add(sp, up, chunk, up + skip, MIN(un - skip, chunk));

	  // copy temporary back
	  mpn_copyi(wp, tp, chunk);
	  *wcy = tcy;

	  return;
	}

      // M >= 2

      if (M % 2)
	{
	  // F(-1) = u_0 - u_1
	  *wcy = -mpn_sub(wp, up, chunk, up + skip, MIN(un - skip, chunk));
	}

      // F(1) = u_0 + u_1
      *scy = mpn_add(sp, up, chunk, up + skip, MIN(un - skip, chunk));

      // recurse into even part
      fwd_impl_1(rp, rcy, up, un, M / 2, 2 * skip, 2 * skip_cy,
		 chunk, scratch);

      // recurse into odd part
      fwd_impl_1(rp + skip, rcy + skip_cy, up + skip, MAX(un, skip) - skip,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
  else if (un <= 3*skip)
    {
      // exactly three input coefficients

      mp_limb_t* tp = scratch;
      scratch += chunk;
      mp_limb_signed_t tcy;

      // recurse into even part, and compute u_0 + u_2
      fwd_impl_2(rp, rcy, tp, &tcy, up, un, M / 2, 2 * skip, 2 * skip_cy,
		 chunk, scratch);

      if (M % 2)
	{
	  // F(-1) = u_0 + u_2 - u_1
	  *wcy = tcy - mpn_sub_n(wp, tp, up + skip, chunk);
	}

      // F(1) = u_0 + u_2 + u_1
      *scy = tcy + mpn_add_n(sp, tp, up + skip, chunk);

      // recurse into odd part
      fwd_impl_1(rp + skip, rcy + skip_cy, up + skip, MAX(un, skip) - skip,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
  else
    {
      // at least four input coefficients

      // recurse into even part, and compute F0(1)
      mp_limb_t* sp2 = scratch;
      scratch += chunk;
      mp_limb_signed_t scy2;

      fwd_impl_2(rp, rcy, sp2, &scy2, up, un, M / 2, 2 * skip, 2 * skip_cy,
		 chunk, scratch);

      // recurse into odd part, and compute F1(1)
      fwd_impl_2(rp + skip, rcy + skip_cy, sp, scy,
		 up + skip, MAX(un, skip) - skip,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      if (M % 2)
	{
	  // F(-1) = F0(1) - F1(1)
	  *wcy = scy2 - *scy - mpn_sub_n(wp, sp2, sp, chunk);
	}

      // F(1) = F0(1) + F1(1)
      *scy += scy2 + mpn_add_n(sp, sp, sp2, chunk);
    }
}


static void
fwd_impl_1(mp_limb_t* rp, mp_limb_signed_t* rcy,
	   mp_limb_t* up, size_t un,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch)
{
  if (M == 0)
    // nothing to do
    return;

  if (un == 0)
    {
      // empty input, so all output coefficients are zero
      for (size_t i = 0; i < M; i++)
      	mpn_zero(rp + i * skip, chunk);
      for (size_t i = 0; i < M; i++)
	rcy[i * skip_cy] = 0;

      return;
    }

  if (M % 2 == 0)
    {
      // M even and M >= 2

      // recurse into odd part
      fwd_impl_1(rp + skip, rcy + skip_cy, up + skip, MAX(un, skip) - skip,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      // recurse into even part
      fwd_impl_1(rp, rcy, up, un, M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      return;
    }

  // M odd

  if (un <= 3*skip)
    {
      // at most three input coefficients

      // sp = where we want to put F(-1)
      mp_limb_t* sp = rp + (M-1) * skip;
      mp_limb_signed_t* scy = rcy + (M-1) * skip_cy;

      if (un <= 2*skip)
	{
	  // at most two input coefficients

	  if (un <= skip)
	    {
	      // exactly one input coefficient

	      // compute F(-1) directly
	      if (sp != up)
		mpn_copyi(sp, up, MIN(un, chunk));
	      if (chunk > un)
		mpn_zero(sp + un, chunk - un);
	      *scy = 0;
	    }
	  else
	    {
	      // exactly two input coefficients

	      // compute F(-1) directly
	      *scy = -mpn_sub(sp, up, chunk, up + skip, MIN(un - skip, chunk));
	    }

	  // recurse into even part
	  fwd_impl_1(rp, rcy, up, un, M / 2, 2 * skip, 2 * skip_cy,
		     chunk, scratch);
	}
      else
	{
	  // exactly three input coefficients

	  // recurse into even part, and compute sum of x^0 and x^2 coefficients
	  fwd_impl_2(rp, rcy, sp, scy, up, un,
		     M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

	  // subtract x^1 coefficient to get F(-1)
	  *scy -= mpn_sub_n(sp, sp, up + skip, chunk);
	}

      // recurse into odd part
      fwd_impl_1(rp + skip, rcy + skip_cy,
     		 up + skip, MAX(un, skip) - skip,
      		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      return;
    }

  // at least four input coefficients

  mp_limb_t* sp = scratch;
  scratch += chunk;
  mp_limb_signed_t scy;

  // recurse into odd part
  fwd_impl_2(rp + skip, rcy + skip_cy, sp, &scy,
	     up + skip, MAX(un, skip) - skip,
	     M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

  // recurse into even part
  fwd_impl_2(rp, rcy, rp + (M-1) * skip, rcy + (M-1) * skip_cy, up, un,
	     M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

  // F(-1) = F0(1) - F1(1)
  rcy[(M-1) * skip_cy] -= scy + mpn_sub_n(rp + (M-1) * skip,
					  rp + (M-1) * skip, sp, chunk);
}



void zz_fermat_mpn_to_transform(zz_fermat_transform_t T,
				mp_limb_t* op, size_t n, int threads)
{
  size_t t = T->params->t;
  size_t M = T->params->M;
  unsigned w = T->params->w;
  unsigned* e = T->params->e;

  // break up t into pieces of size up to ZZ_FERMAT_SLICE

  // number of slices
  size_t slices = (t - 1) / ZZ_FERMAT_SLICE + 1;

  // bound on recursion depth for fwd_impl_*() calls
  unsigned depth = fft62_log2(n / ZZ_FERMAT_SLICE) + 2;

  // cy_buf = "slices" arrays of size M
  mp_limb_signed_t* cy_buf = zz_malloc(M * slices * sizeof(mp_limb_t));

  // main transform
//#pragma omp parallel num_threads(threads)
  {
    mp_limb_t* scratch = (mp_limb_t*)
      zz_malloc(depth * ZZ_FERMAT_SLICE * sizeof(mp_limb_t));

//#pragma omp for schedule(static)
    for (size_t i = 0; i < slices; i++)
      {
	size_t start = i * ZZ_FERMAT_SLICE;
	size_t chunk = MIN(start + ZZ_FERMAT_SLICE, t) - start;
	fwd_impl_1(T->data + start, cy_buf + M * i,
		   op + start, MAX(start, n) - start,
		   M, t, 1, chunk, scratch);
      }

    zz_free(scratch, depth * ZZ_FERMAT_SLICE * sizeof(mp_limb_t));
  }

  // propagate carries
  size_t Msum = 0;
  for (unsigned i = 0; i < w; i++)
    {
      size_t Mi = (size_t) 1 << e[i];
      mp_limb_t* res = T->data + t * Msum;
      mp_limb_t* ptr = res;
      mp_limb_signed_t cy = 0;

      for (size_t j = 0; j < Mi; j++)
	for (size_t k = 0; k < slices; k++)
	  {
	    size_t start = k * ZZ_FERMAT_SLICE;
	    size_t chunk = MIN(start + ZZ_FERMAT_SLICE, t) - start;

	    if (cy >= 0)
	      cy = mpn_add_1(ptr, ptr, chunk, cy);
	    else
	      cy = -mpn_sub_1(ptr, ptr, chunk, -cy);

	    cy += cy_buf[M*k + Msum + j];
	    ptr += chunk;
	  }

      // reduce mod B^(t*M_i) + 1
      if (cy >= 0)
	T->hi[i] = -mpn_sub_1(res, res, t * Mi, cy);
      else
	T->hi[i] = mpn_add_1(res, res, t * Mi, -cy);

      Msum += Mi;
    }

  zz_free(cy_buf, M * slices * sizeof(mp_limb_t));
}


void zz_fermat_mpz_to_transform(zz_fermat_transform_t T, mpz_t op, int threads)
{
  int sign = mpz_sgn(op);

  if (sign == 0)
    zz_fermat_transform_zero(T, threads);
  else
    {
      zz_fermat_mpn_to_transform(T, op->_mp_d, mpz_size(op), threads);
      if (sign == -1)
	zz_fermat_transform_neg(T, T, threads);
    }
}



/*
  inv_impl_* perform inverse fermat transform on polynomials in Z[x].

  Coefficients are represented as a main part ("chunk" limbs), plus a signed
  high limb (weight B^chunk) and a low limb (weight B^(-1)).

  The main parts of each coefficient are spaced by "skip" limbs. The low/high
  limbs are stored in pairs in ucy, spaced by "skip_cy" limbs.

  scratch = scratch space of size ceil(log2(M)) * chunk.

  inv_impl_1(): let F(x) have degree < M. Assumes fermat residues are in up.
  Writes coefficients of F to rp.

  inv_impl_2(): same as inv_impl_1(), but also computes F(1), writes it to sp.

  inv_impl_3(): let F(x) have degree < M + 1. Assumes fermat residues are in up,
  and that M-th coefficient of rp contains F(1). Writes coefficients of F to rp.

  In the comments we write F0(x) and F1(x) for the even and off parts of F,
  i.e. F(x) = F0(x^2) + x F1(x^2).
*/
static void
inv_impl_1(mp_limb_t* rp, mp_limb_t* rcy, mp_limb_t* up, mp_limb_t* ucy,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch);

static void
inv_impl_2(mp_limb_t* rp, mp_limb_t* rcy, mp_limb_t* up, mp_limb_t* ucy,
	   mp_limb_t* sp, mp_limb_t* scy,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch);

static void
inv_impl_3(mp_limb_t* rp, mp_limb_t* rcy, mp_limb_t* up, mp_limb_t* ucy,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch);



static void
inv_coeff_set(mp_limb_t* rp, mp_limb_t* rcy, mp_limb_t* up, mp_limb_t* ucy,
	      size_t chunk)
{
  rcy[0] = ucy[0];
  rcy[1] = ucy[1];
  if (rp != up)
    mpn_copyi(rp, up, chunk);
}


static void
inv_coeff_add(mp_limb_t* rp, mp_limb_t* rcy,
	      mp_limb_t* up1, mp_limb_t* ucy1,
	      mp_limb_t* up2, mp_limb_t* ucy2,
	      size_t chunk)
{
  mp_limb_t cy_lo = mpn_add_1(rcy, ucy1, 1, ucy2[0]);    // low limb
  mp_limb_t cy_hi = mpn_add_n(rp, up1, up2, chunk);      // main chunk
  // high limb
  rcy[1] = ucy1[1] + ucy2[1] + cy_hi + mpn_add_1(rp, rp, chunk, cy_lo);
}


static void
inv_coeff_sub(mp_limb_t* rp, mp_limb_t* rcy,
	      mp_limb_t* up1, mp_limb_t* ucy1,
	      mp_limb_t* up2, mp_limb_t* ucy2,
	      size_t chunk)
{
  mp_limb_t cy_lo = mpn_sub_1(rcy, ucy1, 1, ucy2[0]);    // low limb
  mp_limb_t cy_hi = mpn_sub_n(rp, up1, up2, chunk);      // main chunk
  // high limb
  rcy[1] = ucy1[1] - ucy2[1] - cy_hi - mpn_sub_1(rp, rp, chunk, cy_lo);
}


static void
inv_coeff_rshift(mp_limb_t* rp, mp_limb_t* rcy,
		 mp_limb_t* up, mp_limb_t* ucy,
		 size_t chunk)
{
  rcy[0] = (ucy[0] >> 1) + mpn_rshift(rp, up, chunk, 1);
  rp[chunk - 1] += ucy[1] << 63;
  rcy[1] = (ucy[1] >> 1) + (ucy[1] & ((mp_limb_t) 1 << 63));   // preserve sign
}



// todo: in zz_fermat_transform_to_mpn(), we always operate inplace on
// the cy buffers. So maybe we can remove a few extraneous parameters in
// the following functions?

static void
inv_impl_1(mp_limb_t* rp, mp_limb_t* rcy, mp_limb_t* up, mp_limb_t* ucy,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch)
{
  if (M == 1)
    {
      inv_coeff_set(rp, rcy, up, ucy, chunk);
    }
  else if (M % 2 == 0)
    {
      // recurse into odd part
      inv_impl_1(rp + skip, rcy + skip_cy, up + skip, ucy + skip_cy,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      // recurse into even part
      inv_impl_1(rp, rcy, up, ucy,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
  else if (M == 3)
    {
      // special case to avoid unnecessary copies

      inv_coeff_set(rp + skip, rcy + skip_cy,
		    up + skip, ucy + skip_cy, chunk);

      // F0(1) = F(-1) + F1(1)
      inv_coeff_add(rp + 2 * skip, rcy + 2 * skip_cy,
		    up + 2 * skip, ucy + 2 * skip_cy,
		    up + skip, ucy + skip_cy, chunk);

      // recurse into even part
      inv_impl_3(rp, rcy, up, ucy, 1, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
  else
    {
      // M odd and M >= 5

      mp_limb_t* sp = scratch;
      scratch += chunk;
      mp_limb_t scy[2];

      // recurse into odd part, and compute F1(1)
      inv_impl_2(rp + skip, rcy + skip_cy, up + skip, ucy + skip_cy,
		 sp, scy, M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      // F0(1) = F(-1) + F1(1)
      inv_coeff_add(rp + (M-1) * skip, rcy + (M-1) * skip_cy,
		    up + (M-1) * skip, ucy + (M-1) * skip_cy,
		    sp, scy, chunk);

      // recurse into even part
      inv_impl_3(rp, rcy, up, ucy,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
}


static void
inv_impl_2(mp_limb_t* rp, mp_limb_t* rcy, mp_limb_t* up, mp_limb_t* ucy,
	   mp_limb_t* sp, mp_limb_t* scy,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch)
{
  if (M == 1)
    {
      inv_coeff_set(rp, rcy, up, ucy, chunk);
      inv_coeff_set(sp, scy, up, ucy, chunk);
    }
  else if (M == 2)
    {
      // special case to avoid unnecessary copies

      inv_coeff_set(rp + skip, rcy + skip_cy,
		    up + skip, ucy + skip_cy, chunk);
      inv_coeff_set(rp, rcy, up, ucy, chunk);

      // F(1) = F0(1) + F1(1)
      inv_coeff_add(sp, scy, rp, rcy,
		    rp + skip, rcy + skip_cy, chunk);
    }
  else if (M == 3)
    {
      // special case to avoid unnecessary copies

      inv_coeff_set(rp + skip, rcy + skip_cy,
		    up + skip, ucy + skip_cy, chunk);

      // F0(1) = F(-1) + F1(1)
      inv_coeff_add(rp + 2 * skip, rcy + 2 * skip_cy,
		    up + 2 * skip, ucy + 2 * skip_cy,
		    rp + skip, rcy + skip_cy, chunk);

      // F(1) = F0(1) + F1(1)
      inv_coeff_add(sp, scy, rp + 2 * skip, rcy + 2 * skip_cy,
		    rp + skip, rcy + skip_cy, chunk);

      // recurse into even part
      inv_impl_3(rp, rcy, up, ucy,
		 1, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
  else if (M % 2 == 0)
    {
      // M even and M >= 4

      // recurse into odd part, and compute F1(1)
      inv_impl_2(rp + skip, rcy + skip_cy, up + skip, ucy + skip_cy,
		 sp, scy, M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      mp_limb_t* sp2 = scratch;
      scratch += chunk;
      mp_limb_t scy2[2];

      // recurse into even part, and compute F0(1)
      inv_impl_2(rp, rcy, up, ucy,
		 sp2, scy2, M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      // F(1) = F0(1) + F1(1)
      inv_coeff_add(sp, scy, sp, scy, sp2, scy2, chunk);
    }
  else
    {
      // M odd and M >= 5

      // recurse into odd part, and compute F1(1)
      inv_impl_2(rp + skip, rcy + skip_cy, up + skip, ucy + skip_cy,
		 sp, scy, M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      // F0(1) = F(-1) + F1(1)
      inv_coeff_add(rp + (M-1) * skip, rcy + (M-1) * skip_cy,
		    up + (M-1) * skip, ucy + (M-1) * skip_cy,
		    sp, scy, chunk);

      // F(1) = F0(1) + F1(1)
      inv_coeff_add(sp, scy, sp, scy,
		    rp + (M-1) * skip, rcy + (M-1) * skip_cy, chunk);

      // recurse into even part
      inv_impl_3(rp, rcy, up, ucy,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
}


static void
inv_impl_3(mp_limb_t* rp, mp_limb_t* rcy, mp_limb_t* up, mp_limb_t* ucy,
	   size_t M, size_t skip, size_t skip_cy,
	   size_t chunk, mp_limb_t* scratch)
{
  if (M == 2)
    {
      // special case to avoid unnecessary copies

      inv_coeff_set(rp + skip, rcy + skip_cy,
		    up + skip, ucy + skip_cy, chunk);

      // F0(1) = F(1) - F1(1)
      inv_coeff_sub(rp + 2 * skip, rcy + 2 * skip_cy,
		    rp + 2 * skip, rcy + 2 * skip_cy,
		    rp + skip, rcy + skip_cy, chunk);

      // recurse into even part
      inv_impl_3(rp, rcy, up, ucy,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
  else if (M % 2 == 0)
    {
      // M even and M >= 4

      mp_limb_t* sp = scratch;
      scratch += chunk;
      mp_limb_t scy[2];

      // recurse into odd part, and compute F1(1)
      inv_impl_2(rp + skip, rcy + skip_cy, up + skip, ucy + skip_cy,
		 sp, scy, M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);

      // F0(1) = F(1) - F1(1)
      inv_coeff_sub(rp + M * skip, rcy + M * skip_cy,
		    rp + M * skip, rcy + M * skip_cy,
		    sp, scy, chunk);

      // recurse into even part
      inv_impl_3(rp, rcy, up, ucy,
		 M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
    }
  else
    {
      // M odd and M >= 1

      // F(1) - F(-1)
      inv_coeff_sub(rp + M * skip, rcy + M * skip_cy,
		    rp + M * skip, rcy + M * skip_cy,
		    up + (M-1) * skip, ucy + (M-1) * skip_cy, chunk);

      // (F(1) - F(-1)) / 2
      inv_coeff_rshift(rp + M * skip, rcy + M * skip_cy,
		       rp + M * skip, rcy + M * skip_cy, chunk);

      // (F(1) - F(-1)) / 2 + F(-1)  =  (F(1) + F(-1)) / 2
      inv_coeff_add(rp + (M-1) * skip, rcy + (M-1) * skip_cy,
		    up + (M-1) * skip, ucy + (M-1) * skip_cy,
		    rp + M * skip, rcy + M * skip_cy, chunk);

      // now have F mod x^2 - 1

      if (M > 1)
	{
	  // recurse into odd part
	  inv_impl_3(rp + skip, rcy + skip_cy, up + skip, ucy + skip_cy,
		     M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
	  // recurse into even part
	  inv_impl_3(rp, rcy, up, ucy,
		     M / 2, 2 * skip, 2 * skip_cy, chunk, scratch);
	}
    }
}



void zz_fermat_add_F(mp_limb_t* op, mp_limb_signed_t c,
		     zz_fermat_params_t* params)
{
  size_t t = params->t;
  size_t M = params->M;

  if (c > 0)
    {
      mp_limb_t cy = 0;
      for (size_t i = 0; i < M; i++)
	cy = mpn_add_1(op + i * t, op + i * t, t,
		       cy + (((i & M) == i) ? c : 0));
      op[t * M] += c + cy;
    }
  else if (c < 0)
    {
      mp_limb_t cy = 0;
      for (size_t i = 0; i < M; i++)
	cy = mpn_sub_1(op + i * t, op + i * t, t,
		       cy + (((i & M) == i) ? (-c) : 0));
      op[t * M] += c - cy;
    }
}


// adds c/B mod F to {op, t*M+1}
static void zz_fermat_add_cB(mp_limb_t* op, mp_limb_t c,
			     zz_fermat_params_t* params)
{
  // use the fact that if F = \sum_{i in I} B^i,
  // then -1/B = \sum_{i in I, i != 0} B^{i-1}  mod F

  size_t t = params->t;
  size_t M = params->M;

  mp_limb_t cy = 0;
  for (size_t i = 1; i < M; i++)
    cy = mpn_sub_1(op + i * t - 1, op + i * t - 1, t,
		   cy + (((i & M) == i) ? c : 0));
  mpn_sub_1(op + M * t - 1, op + M * t - 1, 2, cy + c);
}



// todo: the next function doesn't get good parallel speedup for threads > 1.
// Probably either a RAM bandwidth problem, or the threads just take too long to
// get started, but it might be worth investigating further.

void zz_fermat_transform_to_mpn(mp_limb_t* rop, zz_fermat_transform_t T,
				int threads)
{
  size_t t = T->params->t;
  size_t M = T->params->M;
  unsigned w = T->params->w;
  unsigned* e = T->params->e;

  // break up t into pieces of size up to ZZ_FERMAT_SLICE

  // number of slices
  size_t slices = (t - 1) / ZZ_FERMAT_SLICE + 1;

  // cy_buf = "slices" arrays of size 2*M
  mp_limb_t* cy_buf = (mp_limb_t*) zz_malloc(2 * M * slices * sizeof(mp_limb_t));
  mpn_zero(cy_buf, 2 * M * slices);

  // copy high limbs into carries for last slice
  {
    size_t x = 0;
    mp_limb_t* last_slice_cy = cy_buf + 2 * M * (slices - 1);
    for (unsigned i = 0; i < w; i++)
      {
	x += (size_t) 1 << e[i];
	last_slice_cy[2*x - 1] = T->hi[i];
      }
  }

  // main transform
//#pragma omp parallel num_threads(threads)
  {
    mp_limb_t* scratch = (mp_limb_t*)
      zz_malloc(ZZ_FERMAT_MAX_WEIGHT * ZZ_FERMAT_SLICE * sizeof(mp_limb_t));

//#pragma omp for schedule(static)
    for (size_t i = 0; i < slices; i++)
      {
	size_t start = i * ZZ_FERMAT_SLICE;
	size_t chunk = MIN(start + ZZ_FERMAT_SLICE, t) - start;
	inv_impl_1(rop + start, cy_buf + 2 * M * i,
		   T->data + start, cy_buf + 2 * M * i,
		   M, t, 2, chunk, scratch);
      }

    zz_free(scratch, ZZ_FERMAT_MAX_WEIGHT * ZZ_FERMAT_SLICE * sizeof(mp_limb_t));
  }

  // apply carries (low and high limbs)
  mp_limb_signed_t cy = 0;  // carry to apply at current chunk
  mp_limb_t* ptr = rop;     // start of current chunk
  for (size_t j = 0; j < M; j++)
    for (size_t i = 0; i < slices; i++)
      {
	mp_limb_t lo = cy_buf[2*M*i + 2*j];
	mp_limb_signed_t hi = cy_buf[2*M*i + 2*j + 1];

	if (i || j)
	  cy += mpn_add_1(ptr - 1, ptr - 1, 1, lo);

	size_t start = i * ZZ_FERMAT_SLICE;
	size_t chunk = MIN(start + ZZ_FERMAT_SLICE, t) - start;

	if (cy >= 0)
	  cy = mpn_add_1(ptr, ptr, chunk, cy);
	else
	  cy = -mpn_sub_1(ptr, ptr, chunk, -cy);

	cy += hi;
	ptr += chunk;
      }
  rop[t*M] = cy;

  // add in missed low limb, i.e. add cy_buf[0] / B mod F
  zz_fermat_add_cB(rop, cy_buf[0], T->params);

  zz_free(cy_buf, 2 * M * slices * sizeof(mp_limb_t));

  // at this stage the result is {rp,t*M+1}, but still needs to be normalised.

  zz_fermat_add_F(rop, -rop[t*M], T->params);
  assert(rop[t*M] == 0 || rop[t*M] == 1 || (mp_limb_signed_t) rop[t*M] == -1);

  if (rop[t*M] == 1)
    zz_fermat_add_F(rop, -1, T->params);
  assert(rop[t*M] == 0 || (mp_limb_signed_t) rop[t*M] == -1);

  if (rop[t*M])
    zz_fermat_add_F(rop, 1, T->params);
}


void zz_fermat_transform_to_mpz(mpz_t rop, zz_fermat_transform_t T,
				int threads)
{
  size_t n = T->params->t * T->params->M + 1;
  mpz_realloc(rop, n);
  zz_fermat_transform_to_mpn(rop->_mp_d, T, threads);
  // normalise
  while (n > 0 && rop->_mp_d[n - 1] == 0)
    n--;
  rop->_mp_size = n;
}



void zz_fermat_transform_zero(zz_fermat_transform_t T, int threads)
{
  mpn_zero(T->data, T->params->M * T->params->t);
  for (unsigned i = 0; i < T->params->w; i++)
    T->hi[i] = 0;
}


void
zz_fermat_transform_sub(zz_fermat_transform_t rop, zz_fermat_transform_t op1,
			zz_fermat_transform_t op2, int threads)
{
  assert(op1->params == rop->params);
  assert(op2->params == rop->params);

  size_t t = rop->params->t;
  unsigned w = rop->params->w;
  unsigned* e = rop->params->e;

  mp_limb_t* op1_ptr = op1->data;
  mp_limb_t* op2_ptr = op2->data;
  mp_limb_t* rop_ptr = rop->data;

  for (unsigned i = 0; i < w; i++)
    {
      // subtract modulo B^{M_i*t} + 1
      size_t len = t << e[i];     // M_i * t

      mp_limb_signed_t hi = op1->hi[i] - op2->hi[i] -
	mpn_sub_n(rop_ptr, op1_ptr, op2_ptr, len);

      if (hi == 2)
	hi = 1 - mpn_sub_1(rop_ptr, rop_ptr, len, 1);
      else if (hi == (mp_limb_signed_t) (-2))
	hi = mpn_add_1(rop_ptr, rop_ptr, len, 1) - 1;

      rop->hi[i] = hi;

      op1_ptr += len;
      op2_ptr += len;
      rop_ptr += len;
    }
}


void
zz_fermat_transform_neg(zz_fermat_transform_t rop, zz_fermat_transform_t op,
			int threads)
{
  assert(op->params == rop->params);

  size_t t = rop->params->t;
  unsigned w = rop->params->w;
  unsigned* e = rop->params->e;

  mp_limb_t* op_ptr = op->data;
  mp_limb_t* rop_ptr = rop->data;

  for (unsigned i = 0; i < w; i++)
    {
      // negate modulo B^{M_i*t} + 1
      size_t len = t << e[i];     // M_i * t

      mp_limb_signed_t hi = -op->hi[i] - mpn_neg(rop_ptr, op_ptr, len);

      if (hi == (mp_limb_signed_t) (-2))
	hi = mpn_add_1(rop_ptr, rop_ptr, len, 1) - 1;

      rop->hi[i] = hi;

      op_ptr += len;
      rop_ptr += len;
    }
}



void zz_fermat_poly_init(zz_fermat_poly_t P, zz_fermat_params_t* params)
{
  P->params = params;
  P->data[0][0] = NULL;
}


void zz_fermat_poly_clear(zz_fermat_poly_t P)
{
  zz_fermat_poly_dealloc(P);
}


void zz_fermat_poly_alloc(zz_fermat_poly_t P)
{
  if (P->data[0][0] == NULL)
    {
      size_t M = P->params->M;
      size_t K = P->params->K;
      unsigned num_primes = P->params->num_primes;
      unsigned w = P->params->w;
      unsigned* e = P->params->e;

      // allocate in one big chunk
      uint64_t* ptr = (uint64_t*) zz_malloc(num_primes * M * K * sizeof(uint64_t));

      // set up pointers
      for (unsigned j = 0; j < num_primes; j++)
	for (unsigned i = 0; i < w; i++)
	  {
	    P->data[i][j] = ptr;
	    ptr += K << e[i];
	  }
    }
}


void zz_fermat_poly_dealloc(zz_fermat_poly_t P)
{
  if (P->data[0][0] != NULL)
    {
      zz_free(P->data[0][0], P->params->num_primes * P->params->M * P->params->K * sizeof(uint64_t));
      P->data[0][0] = NULL;
    }
}


void zz_fermat_transform_to_poly(zz_fermat_poly_t P, zz_fermat_transform_t T,
				 int sign, int lgS, int scale, int threads)
{
  assert(P->params == T->params);

  zz_moduli_t* moduli = P->params->moduli;
  unsigned num_primes = P->params->num_primes;
  unsigned w = P->params->w;
  unsigned lgK = P->params->lgK;
  size_t K = P->params->K;
  unsigned* e = P->params->e;
  unsigned r = P->params->r;
  size_t t = P->params->t;

  zz_fermat_poly_alloc(P);

  mp_limb_t* src = T->data;
  mp_limb_signed_t* hi = T->hi;

  for (unsigned i = 0; i < w; i++)
    {
      uint64_t** dst = P->data[i];
      int lgC = lgS - (scale ? (lgK + e[i]) : 0);

      zz_split_reduce(dst, K << e[i], sign, lgC, moduli, num_primes,
		      src, t << e[i], r, threads);

      // contribution from overflow limb
      // (subtract it from the constant term, modulo each prime separately)
      if (hi[i])
	for (unsigned j = 0; j < num_primes; j++)
	  {
	    uint64_t p = moduli->p[j];
	    uint64_t pinv = moduli->pinv[j];
	    uint64_t c = mod62_2exp(lgC, p, pinv);
	    c = (sign == (int) hi[i]) ? c : (p - c);
	    dst[j][0] = mod63_sub(dst[j][0], c, 2*p);
	  }

      src += t << e[i];
    }
}


void zz_fermat_poly_to_transform(zz_fermat_transform_t T, zz_fermat_poly_t P,
				 int threads)
{
  assert(P->params == T->params);

  zz_moduli_t* moduli = P->params->moduli;
  unsigned num_primes = P->params->num_primes;
  unsigned w = P->params->w;
  size_t K = P->params->K;
  unsigned* e = P->params->e;
  unsigned r = P->params->r;
  size_t t = P->params->t;

  mp_limb_t* dst = T->data;
  mp_limb_signed_t* hi = T->hi;

  for (unsigned i = 0; i < w; i++)
    {
      // number of limbs needed for evaluation
      // (this will leak into the next residue, but it doesn't matter, because
      // we work left-to-right, reducing as we go, and the last one is covered
      // by ZZ_FERMAT_PADDING)
      size_t len = (((K << e[i]) - 1) * r + 62 * num_primes + 2) / 64 + 1;

      zz_crt_recompose(dst, len, r, P->data[i], K << e[i],
		       moduli, num_primes, threads);

      // reduce modulo B^(t*M_i) + 1
      size_t d = t << e[i];
      mp_limb_t borrow = dst[len - 1] >> 63;
      for (size_t j = len - 1; j >= d; j--)
	{
	  if (borrow)
	    borrow = 1 - mpn_add_1(dst + j - d, dst + j - d, d, -dst[j] - 1);
	  else
	    borrow = mpn_sub_1(dst + j - d, dst + j - d, d, dst[j]);
	}
      hi[i] = -borrow;

      dst += t << e[i];
    }
}



void zz_fermat_poly_fft(zz_fermat_poly_t rop, zz_fermat_poly_t op,
			int threads)
{
  assert(op->params == rop->params);

  zz_fermat_poly_alloc(rop);

  unsigned num_primes = rop->params->num_primes;
  unsigned w = rop->params->w;
  unsigned* e = rop->params->e;
  unsigned lgK = rop->params->lgK;

  unsigned teams = zz_gcd(threads, num_primes);
  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned j = 0; j < num_primes; j++)
    {
      fft62_mod_t* mod = &rop->params->moduli->fft62_mod[j];

      for (unsigned i = 0; i < w; i++)
	{
	  unsigned lgN = e[i] + lgK;
	  size_t N = (size_t) 1 << lgN;
	  uint64_t* dst = rop->data[i][j];
	  uint64_t* src = op->data[i][j];
	  fft62_fft_twisted(dst, N, src, N, lgN, 0, lgN + 1, mod, threads2);
	}
    }
}


void zz_fermat_poly_ifft(zz_fermat_poly_t rop, zz_fermat_poly_t op,
			 int threads)
{
  assert(op->params == rop->params);

  zz_fermat_poly_alloc(rop);

  unsigned num_primes = rop->params->num_primes;
  unsigned w = rop->params->w;
  unsigned* e = rop->params->e;
  unsigned lgK = rop->params->lgK;

  unsigned teams = zz_gcd(threads, num_primes);
  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned j = 0; j < num_primes; j++)
    {
      fft62_mod_t* mod = &rop->params->moduli->fft62_mod[j];

      for (unsigned i = 0; i < w; i++)
	{
	  uint64_t* src = op->data[i][j];
	  uint64_t* dst = rop->data[i][j];
	  unsigned lgN = e[i] + lgK;
	  size_t N = (size_t) 1 << lgN;

	  fft62_ifft_twisted(dst, N, src, lgN, 0, lgN + 1, mod, threads2);
	}
    }
}


void zz_fermat_poly_mul(zz_fermat_poly_t rop, zz_fermat_poly_t op1,
			zz_fermat_poly_t op2, int use_pinvb, int threads)
{
  assert(op1->params == rop->params);
  assert(op2->params == rop->params);

  zz_fermat_poly_alloc(rop);

  unsigned num_primes = rop->params->num_primes;
  unsigned w = rop->params->w;
  unsigned* e = rop->params->e;
  unsigned lgK = rop->params->lgK;

//  unsigned teams = zz_gcd(threads, num_primes);
//  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned j = 0; j < num_primes; j++)
    {
      fft62_mod_t* mod = &rop->params->moduli->fft62_mod[j];

      for (unsigned i = 0; i < w; i++)
	{
	  uint64_t p = mod->p;
	  uint64_t* src1 = op1->data[i][j];
	  uint64_t* src2 = op2->data[i][j];
	  uint64_t* dst = rop->data[i][j];
	  unsigned lgN = e[i] + lgK;
	  size_t N = (size_t) 1 << lgN;

	  if (use_pinvb)
	    {
	      uint64_t pinvb = mod->pinvb;
//#pragma omp parallel for num_threads(threads2) schedule(static)
	      for (size_t h = 0; h < N; h++)
		dst[h] = mod62_mul_pinvb_lazy(src1[h], src2[h], p, pinvb);
	    }
	  else
	    {
	      uint64_t pinv = mod->pinv;
//#pragma omp parallel for num_threads(threads2) schedule(static)
	      for (size_t h = 0; h < N; h++)
		dst[h] = mod62_mul_pinv_lazy2(src1[h], src2[h], p, pinv);
	    }
	}
    }
}


void zz_fermat_poly_add(zz_fermat_poly_t rop, zz_fermat_poly_t op1,
			zz_fermat_poly_t op2, int threads)
{
  assert(op1->params == rop->params);
  assert(op2->params == rop->params);

  zz_fermat_poly_alloc(rop);

  unsigned num_primes = rop->params->num_primes;
  unsigned w = rop->params->w;
  unsigned* e = rop->params->e;
  unsigned lgK = rop->params->lgK;

//  unsigned teams = zz_gcd(threads, num_primes);
//  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned j = 0; j < num_primes; j++)
    {
      fft62_mod_t* mod = &rop->params->moduli->fft62_mod[j];
      uint64_t p = mod->p;

      for (unsigned i = 0; i < w; i++)
	{
	  uint64_t* src1 = op1->data[i][j];
	  uint64_t* src2 = op2->data[i][j];
	  uint64_t* dst = rop->data[i][j];
	  unsigned lgN = e[i] + lgK;
	  size_t N = (size_t) 1 << lgN;

//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = 0; h < N; h++)
	    dst[h] = mod63_add(src1[h], src2[h], 2*p);
	}
    }
}



void zz_fermat_poly_sub(zz_fermat_poly_t rop, zz_fermat_poly_t op1,
			zz_fermat_poly_t op2, int threads)
{
  assert(op1->params == rop->params);
  assert(op2->params == rop->params);

  zz_fermat_poly_alloc(rop);

  unsigned num_primes = rop->params->num_primes;
  unsigned w = rop->params->w;
  unsigned* e = rop->params->e;
  unsigned lgK = rop->params->lgK;

//  unsigned teams = zz_gcd(threads, num_primes);
//  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned j = 0; j < num_primes; j++)
    {
      fft62_mod_t* mod = &rop->params->moduli->fft62_mod[j];
      uint64_t p = mod->p;

      for (unsigned i = 0; i < w; i++)
	{
	  uint64_t* src1 = op1->data[i][j];
	  uint64_t* src2 = op2->data[i][j];
	  uint64_t* dst = rop->data[i][j];
	  unsigned lgN = e[i] + lgK;
	  size_t N = (size_t) 1 << lgN;

//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = 0; h < N; h++)
	    dst[h] = mod63_sub(src1[h], src2[h], 2*p);
	}
    }
}
