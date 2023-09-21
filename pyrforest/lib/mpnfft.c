/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <assert.h>
#include <math.h>
#include "mpnfft.h"
#include "split_reduce.h"
#include "crt_recompose.h"
#include "zzmisc.h"
#include "fft62/fft62.h"
#include "fft62/mod62.h"
#include "zzmem.h"


void zz_mpnfft_params_init(zz_mpnfft_params_t* params, size_t bits,
			   size_t terms, unsigned num_primes,
			   zz_moduli_t* moduli)
{
  params->moduli = moduli;
  params->num_primes = num_primes;
  params->terms = terms;

  // we need that if |x1*y1| + ... + |xterms*yterms| fits into "bits" bits, then
  // "points" is big enough to accommodate product of polynomials representing
  // xi and yi.

  // so we need: if terms * |xi| * |yi| < 2^bits, then
  //   points >= ceil(bits(|xi|)/r) + ceil(bits(|yi|)/r) - 1.

  // after some algebra, we find it suffices to get points >= bits / r + 1.

  // the other constraint is that we need the coefficients of the product
  // polynomial to fit into the product of the FFT primes, i.e.
  //    terms * ceil(points / 2) * 2^(2r) < (product of FFT primes) / 2.
  
  // (here ceil(points / 2) is an upper bound for the length of the shorter of
  // the polynomials representing xi and yi)

  // bound on polynomial coefficients
  // todo: this will overflow IEEE double after about 16 primes or so
  double bound = 0.4999;
  for (unsigned i = 0; i < num_primes; i++)
    bound *= (double) moduli->p[i];
  double bound2 = bound / terms;

  // initial (over)estimate of r
  unsigned r = 62 * num_primes / 2;

  // this loop generates a decreasing sequence of overestimates of r
  // until we hit one that works
  while (1)
    {
      size_t points = (bits - 1) / r + 2;    // ceil(bits / r + 1)

      // find largest r' such that 2^(2r') < bound2 / ceil(points / 2)
      double thing = bound2 / ((points + 1) / 2);
      unsigned rr;
      frexp(thing, (int*) &rr);
      rr = (rr - 1) / 2;

      // if r' matches r, we're done
      if (rr == r)
	{
	  // done
	  unsigned lgN = fft62_log2(points);
	  params->r = r;
	  params->lgN = lgN; 
	  params->N = (size_t) 1 << lgN;
	  params->points = fft62_next_size(points, lgN);
	  return;
	}
      // otherwise decrease r and try again
      r = rr;
    }
}


void zz_mpnfft_poly_init(zz_mpnfft_poly_t P, zz_mpnfft_params_t* params)
{
  P->params = params;
  P->size = 0;
  P->data[0] = NULL;
}


void zz_mpnfft_poly_clear(zz_mpnfft_poly_t P)
{
  zz_mpnfft_poly_dealloc(P);
}


void zz_mpnfft_poly_alloc(zz_mpnfft_poly_t P)
{
  if (P->data[0] == NULL)
    {
      size_t points = P->params->points;
      unsigned num_primes = P->params->num_primes;
      P->data[0] = (uint64_t*) zz_malloc(num_primes * points * sizeof(uint64_t));
      for (unsigned i = 1; i < num_primes; i++)
	P->data[i] = P->data[i - 1] + points;
    }
}


void zz_mpnfft_poly_dealloc(zz_mpnfft_poly_t P)
{
  if (P->data[0] != NULL)
    {
      zz_free(P->data[0], P->params->num_primes * P->params->points * sizeof(uint64_t));
      P->data[0] = NULL;
    }
  P->size = 0;
}


void zz_mpnfft_mpn_to_poly(zz_mpnfft_poly_t P, mp_limb_t* up, size_t un,
			   int sign, int lgS, int scale, int threads)
{
  if (sign == 0 || un == 0)
    {
      P->size = 0;
      return;
    }

  zz_mpnfft_poly_alloc(P);

  unsigned lgN = P->params->lgN;
  zz_moduli_t* moduli = P->params->moduli;
  unsigned num_primes = P->params->num_primes;
  size_t points = P->params->points;
  unsigned r = P->params->r;

  // estimate number of coefficients needed
  size_t size = (64 * un + r - 1) / r;
  if (size > points)
    size = points;
  size = fft62_next_size(size, P->params->lgN);
  P->size = size;
  
  zz_split_reduce(P->data, size, sign, lgS - (scale ? lgN : 0),
		  moduli, num_primes, up, un, r, threads);
}


void zz_mpnfft_poly_to_mpn(mp_limb_t* rp, size_t rn, zz_mpnfft_poly_t P,
			   int threads)
{
  zz_crt_recompose(rp, rn, P->params->r, P->data, P->size, P->params->moduli,
		   P->params->num_primes, threads);
}


void zz_mpnfft_poly_set(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op, int threads)
{
  assert(rop->params == op->params);

  if (rop == op)
    return;

  if (op->size == 0)
    {
      rop->size = 0;
      return;
    }

  unsigned num_primes = rop->params->num_primes;
  size_t size = op->size;

  zz_mpnfft_poly_alloc(rop);

//  unsigned teams = zz_gcd(threads, num_primes);
//  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned i = 0; i < num_primes; i++)
    {
      uint64_t* src = op->data[i];
      uint64_t* dst = rop->data[i];

//#pragma omp parallel for num_threads(threads2) schedule(static)
      for (size_t h = 0; h < size; h++)
	dst[h] = src[h];
    }

  rop->size = size;
}


void zz_mpnfft_poly_neg(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op, int threads)
{
  assert(rop->params == op->params);

  if (op->size == 0)
    {
      rop->size = 0;
      return;
    }

  zz_moduli_t* moduli = rop->params->moduli;
  unsigned num_primes = rop->params->num_primes;
  size_t size = op->size;

  zz_mpnfft_poly_alloc(rop);

  // unsigned teams = zz_gcd(threads, num_primes);
  // unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned i = 0; i < num_primes; i++)
    {
      fft62_mod_t* mod = &moduli->fft62_mod[i];
      uint64_t p = mod->p;
      uint64_t* src = op->data[i];
      uint64_t* dst = rop->data[i];

//#pragma omp parallel for num_threads(threads2) schedule(static)
      for (size_t h = 0; h < size; h++)
	dst[h] = mod63_sub(0, src[h], 2*p);
    }

  rop->size = size;
}


void zz_mpnfft_poly_add(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op1,
			zz_mpnfft_poly_t op2, int threads)
{
  assert(rop->params == op1->params);
  assert(rop->params == op2->params);

  if (op1->size == 0)
    {
      zz_mpnfft_poly_set(rop, op2, threads);
      return;
    }

  if (op2->size == 0)
    {
      zz_mpnfft_poly_set(rop, op1, threads);
      return;
    }

  size_t size1 = op1->size;
  size_t size2 = op2->size;
  size_t size_min = MIN(size1, size2);
  size_t size_max = MAX(size1, size2);

  zz_moduli_t* moduli = rop->params->moduli;
  unsigned num_primes = rop->params->num_primes;

  zz_mpnfft_poly_alloc(rop);

//  unsigned teams = zz_gcd(threads, num_primes);
//  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned i = 0; i < num_primes; i++)
    {
      fft62_mod_t* mod = &moduli->fft62_mod[i];
      uint64_t p = mod->p;

      uint64_t* src1 = op1->data[i];
      uint64_t* src2 = op2->data[i];
      uint64_t* dst = rop->data[i];

//#pragma omp parallel for num_threads(threads2) schedule(static)
      for (size_t h = 0; h < size_min; h++)
	dst[h] = mod63_add(src1[h], src2[h], 2*p);

      if (size1 <= size2)
	{
//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = size1; h < size2; h++)
	    dst[h] = src2[h];
	}
      else
	{
//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = size2; h < size1; h++)
	    dst[h] = src1[h];
	}
    }

  rop->size = size_max;
}


void zz_mpnfft_poly_sub(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op1,
			zz_mpnfft_poly_t op2, int threads)
{
  assert(rop->params == op1->params);
  assert(rop->params == op2->params);

  if (op1->size == 0)
    {
      zz_mpnfft_poly_neg(rop, op2, threads);
      return;
    }

  if (op2->size == 0)
    {
      zz_mpnfft_poly_set(rop, op1, threads);
      return;
    }

  size_t size1 = op1->size;
  size_t size2 = op2->size;
  size_t size_min = MIN(size1, size2);
  size_t size_max = MAX(size1, size2);

  zz_moduli_t* moduli = rop->params->moduli;
  unsigned num_primes = rop->params->num_primes;

  zz_mpnfft_poly_alloc(rop);

//  unsigned teams = zz_gcd(threads, num_primes);
//  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned i = 0; i < num_primes; i++)
    {
      fft62_mod_t* mod = &moduli->fft62_mod[i];
      uint64_t p = mod->p;

      uint64_t* src1 = op1->data[i];
      uint64_t* src2 = op2->data[i];
      uint64_t* dst = rop->data[i];

//#pragma omp parallel for num_threads(threads2) schedule(static)
      for (size_t h = 0; h < size_min; h++)
	dst[h] = mod63_sub(src1[h], src2[h], 2*p);

      if (size1 <= size2)
	{
//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = size1; h < size2; h++)
	    dst[h] = mod63_sub(0, src2[h], 2*p);
	}
      else
	{
//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = size2; h < size1; h++)
	    dst[h] = src1[h];
	}
    }

  rop->size = size_max;
}


void zz_mpnfft_poly_fft(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op,
			int threads)
{
  assert(rop->params == op->params);

  if (op->size == 0)
    {
      rop->size = 0;
      return;
    }

  zz_mpnfft_poly_alloc(rop);

  zz_moduli_t* moduli = rop->params->moduli;
  unsigned num_primes = rop->params->num_primes;
  unsigned lgN = rop->params->lgN;
  size_t points = rop->params->points;

  unsigned teams = zz_gcd(threads, num_primes);
  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned i = 0; i < num_primes; i++)
    {
      fft62_fft(rop->data[i], points, op->data[i], op->size, lgN,
		&moduli->fft62_mod[i], threads2);
    }

  rop->size = points;
}


void zz_mpnfft_poly_ifft(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op,
			 int scale, int threads)
{
  assert(rop->params == op->params);

  if (op->size == 0)
    {
      rop->size = 0;
      return;
    }

  zz_mpnfft_poly_alloc(rop);

  zz_moduli_t* moduli = rop->params->moduli;
  unsigned num_primes = rop->params->num_primes;
  unsigned lgN = rop->params->lgN;
  size_t points = rop->params->points;

  assert(op->size == points);

  unsigned teams = zz_gcd(threads, num_primes);
  unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned i = 0; i < num_primes; i++)
    {
      fft62_mod_t* mod = &moduli->fft62_mod[i];
      uint64_t* src = op->data[i];
      uint64_t* dst = rop->data[i];

      fft62_ifft(dst, points, src, lgN, mod, threads2);

      if (scale)
	{
	  // divide by 2^lgN
	  uint64_t p = mod->p;
	  uint64_t pinv = mod->pinv;
	  uint64_t u = mod62_2exp(-lgN, p, pinv);
	  uint64_t upinv = mod62_ypinv(u, p, pinv);
//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = 0; h < points; h++)
	    dst[h] = mod62_mul_ypinv_lazy(dst[h], u, upinv, p);
	}
    }

  rop->size = points;
}


void zz_mpnfft_poly_mul(zz_mpnfft_poly_t rop, zz_mpnfft_poly_t op1,
			zz_mpnfft_poly_t op2, int use_pinvb, int threads)
{
  assert(op1->params == rop->params);
  assert(op2->params == rop->params);

  if (op1->size == 0 || op2->size == 0)
    {
      rop->size = 0;
      return;
    }

  zz_moduli_t* moduli = rop->params->moduli;
  unsigned num_primes = rop->params->num_primes;
  size_t points = rop->params->points;

  assert(op1->size == points);
  assert(op2->size == points);

  zz_mpnfft_poly_alloc(rop);

  // unsigned teams = zz_gcd(threads, num_primes);
  // unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned i = 0; i < num_primes; i++)
    {
      fft62_mod_t* mod = &moduli->fft62_mod[i];
      uint64_t* src1 = op1->data[i];
      uint64_t* src2 = op2->data[i];
      uint64_t* dst = rop->data[i];

      uint64_t p = mod->p;
      if (use_pinvb)
	{
	  uint64_t pinvb = mod->pinvb;
//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = 0; h < points; h++)
	    dst[h] = mod62_mul_pinvb_lazy(src1[h], src2[h], p, pinvb);
	}
      else
	{
	  uint64_t pinv = mod->pinv;
//#pragma omp parallel for num_threads(threads2) schedule(static)
	  for (size_t h = 0; h < points; h++)
	    dst[h] = mod62_mul_pinv_lazy2(src1[h], src2[h], p, pinv);
	}
    }

  rop->size = points;
}



// helper function for zz_mpnfft_poly_matrix_mul()
static void chunk_copy_in(uint64_t* dst, zz_mpnfft_poly_t op, unsigned g,
			  size_t h0, size_t len, uint64_t p)
{
  if (op->size == 0)
    {
      for (size_t h = 0; h < len; h++)
	dst[h] = 0;
    }
  else
    {
      uint64_t* src = op->data[g] + h0;
      for (size_t h = 0; h < len; h++)
	dst[h] = mod63_reduce2(src[h], p);
    }
}


// helper function for zz_mpnfft_poly_matrix_mul()
static void chunk_reduce(uint64_t* dest, size_t len, uint64_t p, uint64_t pinv,
			 uint64_t Bmodp)
{
  for (size_t h = 0; h < len; h++)
    {
      uint64_t u1, u0;
      MUL_WIDE(u1, u0, dest[2*h+1], Bmodp);
      ADD_WIDE(u1, u0, 0, dest[2*h], u1, u0);
      dest[2*h] = mod62_reduce_pinv_lazy4(u1, u0, p, pinv);
      dest[2*h+1] = 0;
    }
}


// helper function for zz_mpnfft_poly_matrix_mul()
static void chunk_reduce_out(uint64_t* dest, uint64_t* src, size_t len,
			     uint64_t p, uint64_t pinv, uint64_t Bmodp)
{
  for (size_t h = 0; h < len; h++)
    {
      uint64_t u1, u0;
      MUL_WIDE(u1, u0, src[2*h+1], Bmodp);
      ADD_WIDE(u1, u0, 0, src[2*h], u1, u0);
      dest[h] = mod62_reduce_pinv_lazy2(u1, u0, p, pinv);
    }
}


// helper function for mpzfft_matrix_mul()
static void chunk_mul(uint64_t* dest, uint64_t* src1, uint64_t* src2,
		      size_t len)
{
  for (size_t h = 0; h < len; h++)
    MUL_WIDE(dest[2*h+1], dest[2*h], src1[h], src2[h]);   // in [0, p^2)
}


// helper function for mpzfft_matrix_mul()
static void chunk_addmul(uint64_t* dest, uint64_t* src1, uint64_t* src2,
			 size_t len)
{
  for (size_t h = 0; h < len; h++)
    {
      uint64_t u1, u0;
      MUL_WIDE(u1, u0, src1[h], src2[h]);   // in [0, p^2)
      ADD_WIDE(dest[2*h+1], dest[2*h], u1, u0, dest[2*h+1], dest[2*h]);
    }
}



void zz_mpnfft_poly_matrix_mul(zz_mpnfft_poly_t* rop,
			       zz_mpnfft_poly_t* op1,
			       zz_mpnfft_poly_t* op2,
			       unsigned dim1, unsigned dim2, unsigned dim3,
			       int threads)
{
  zz_moduli_t* moduli = (*rop)->params->moduli;
  unsigned num_primes = (*rop)->params->num_primes;
  size_t points = (*rop)->params->points;

  for (unsigned i = 0; i < dim1 * dim3; i++)
    zz_mpnfft_poly_alloc(rop[i]);

  // work in chunks, try to fit everything in L1
  size_t chunk = 4096 / (dim1 * dim2 + dim2 * dim3 + dim1 * dim3);
  if (chunk < 16)
    chunk = 16;   // but also don't want inner loop too short

  // unsigned teams = zz_gcd(threads, num_primes);
  // unsigned threads2 = threads / teams;

//#pragma omp parallel for num_threads(teams) schedule(static)
  for (unsigned g = 0; g < num_primes; g++)
    {
      fft62_mod_t* mod = &moduli->fft62_mod[g];
      uint64_t p = mod->p;
      uint64_t pinv = mod->pinv;
      uint64_t Bmodp = mod->B;

//#pragma omp parallel num_threads(threads2)
      {
	uint64_t* src1_temp = zz_malloc(dim1 * dim2 * chunk * sizeof(uint64_t));
	uint64_t* src2_temp = zz_malloc(dim2 * dim3 * chunk * sizeof(uint64_t));
	uint64_t* dest_temp = zz_malloc(2 * dim1 * dim3 * chunk *
				     sizeof(uint64_t));

//#pragma omp for schedule(static)
	for (size_t h0 = 0; h0 < points; h0 += chunk)
	  {
	    // this iteration handles coefficients in h0 <= h < h1
	    size_t h1 = MIN(h0 + chunk, points);
	    size_t len = h1 - h0;

	    // copy input matrices, reduce from [0, 2p) to [0, p)
	    for (unsigned i = 0; i < dim1 * dim2; i++)
	      chunk_copy_in(src1_temp + i * chunk, op1[i], g, h0, len, p);
	    for (unsigned i = 0; i < dim2 * dim3; i++)
	      chunk_copy_in(src2_temp + i * chunk, op2[i], g, h0, len, p);

	    // main matrix multiplication loop
	    for (unsigned i = 0; i < dim1; i++)
	      for (unsigned k = 0; k < dim3; k++)
		{
		  uint64_t* dest = dest_temp + 2 * (i * dim3 + k) * chunk;
		  uint64_t* src1 = src1_temp + i * dim2 * chunk;
		  uint64_t* src2 = src2_temp + k * chunk;

		  chunk_mul(dest, src1, src2, len);

		  for (unsigned j = 1; j < dim2; j++)
		    {
		      src1 += chunk;
		      src2 += dim3 * chunk;

		      // each iteration adds at most p^2 to each entry.
		      // We need to reduce after every 16 iterations, since
		      //   16p^2 + 4p < 2^128
		      if (j % 16 == 0)
			chunk_reduce(dest, len, p, pinv, Bmodp);

		      chunk_addmul(dest, src1, src2, len);
		    }
		}

	    // copy output to rop, and reduce to [0, 2p)
	    for (unsigned i = 0; i < dim1; i++)
	      for (unsigned k = 0; k < dim3; k++)
		{
		  uint64_t* dest = dest_temp + 2 * (i * dim3 + k) * chunk;
		  uint64_t* out = rop[i * dim3 + k]->data[g] + h0;
		  chunk_reduce_out(out, dest, len, p, pinv, Bmodp);
		}
	  }

	zz_free(src1_temp, dim1 * dim2 * chunk * sizeof(uint64_t));
	zz_free(src2_temp, dim2 * dim3 * chunk * sizeof(uint64_t));
	zz_free(dest_temp, 2 * dim1 * dim3 * chunk * sizeof(uint64_t));
      }      
    }

  for (unsigned i = 0; i < dim1 * dim3; i++)
    rop[i]->size = points;
}
