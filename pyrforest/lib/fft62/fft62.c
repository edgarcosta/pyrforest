/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <assert.h>
#include <stddef.h>
#include "mod62.h"
#include "fft62.h"
#include "../zzmem.h"


#ifdef TEST

unsigned FFT62_ARRAY_THRESHOLD = FFT62_ARRAY_THRESHOLD_DEFAULT;
unsigned FFT62_CACHE_THRESHOLD = FFT62_CACHE_THRESHOLD_DEFAULT;

#endif


#define MIN(xxx, yyy) (((xxx) < (yyy)) ? (xxx) : (yyy))
#define MAX(xxx, yyy) (((xxx) > (yyy)) ? (xxx) : (yyy))



unsigned fft62_log2(size_t n)
{
  if (n == 0)
    return 0;
  n--;
  unsigned bits = 0;
  while (n >= 16)
    bits += 4, n >>= 4;
  return bits + ((0x4444444433332210ULL >> (4 * n)) & 15);
}



// finds primitive (2^FFT62_MAX_LGN)-th root mod p
static uint64_t fft62_primitive_root(uint64_t p, uint64_t pinv)
{
  uint64_t N = (uint64_t) 1 << FFT62_MAX_LGN;

  for (uint64_t u = 2; ; u++)
    {
      uint64_t w = mod62_pow_pinv(u, (p - 1) / N, p, pinv);
      if (mod62_pow_pinv(w, N / 2, p, pinv) != 1)
	return w;
    }
}



void fft62_mod_init(fft62_mod_t* mod, uint64_t p)
{
  assert(mod62_valid(p));
  assert(p % ((uint64_t) 1 << FFT62_MAX_LGN) == 1);

  mod->p = p;

  uint64_t pinv = mod->pinv = mod62_pinv(p);
  mod->pinvb = mod62_pinvb(p);

  mod->B = mod62_mul_pinv((uint64_t) 1 << 32, (uint64_t) 1 << 32, p, pinv);

  uint64_t g = mod->g = fft62_primitive_root(p, pinv);
  uint64_t ginv = mod->ginv = mod62_inv(g, p);

  // initialise rtab, rinvtab

  mod->rtab[FFT62_MAX_LGN] = g;
  mod->rinvtab[FFT62_MAX_LGN] = ginv;
  for (int lgN = FFT62_MAX_LGN - 1; lgN >= 0; lgN--)
    {
      mod->rtab[lgN] =
	mod62_mul_pinv(mod->rtab[lgN + 1], mod->rtab[lgN + 1], p, pinv);
      mod->rinvtab[lgN] =
	mod62_mul_pinv(mod->rinvtab[lgN + 1], mod->rinvtab[lgN + 1], p, pinv);
    }

  // initialise root tables

  for (unsigned lgN = 0; lgN <= FFT62_MAX_LGN; lgN++)
    {
      size_t N = (uint64_t) 1 << lgN;

      if (lgN >= 1 && lgN <= FFT62_ARRAY_THRESHOLD)
	{
	  uint64_t w = mod->rtab[lgN];
	  uint64_t wpinv = mod62_ypinv(w, p, pinv);
	  uint64_t winv = mod->rinvtab[lgN];
	  uint64_t winvpinv = mod62_ypinv(winv, p, pinv);

	  uint64_t* wtab = zz_malloc(N / 2 * sizeof(uint64_t));
	  uint64_t* wpinvtab = zz_malloc(N / 2 * sizeof(uint64_t));
	  uint64_t* winvtab = zz_malloc(N / 2 * sizeof(uint64_t));
	  uint64_t* winvpinvtab = zz_malloc(N / 2 * sizeof(uint64_t));

	  mod->wtab[lgN] = wtab;
	  mod->wpinvtab[lgN] = wpinvtab;
	  mod->winvtab[lgN] = winvtab;
	  mod->winvpinvtab[lgN] = winvpinvtab;

	  wtab[0] = 1;
	  winvtab[0] = 1;
	  for (size_t i = 1; i < N / 2; i++)
	    {
	      wtab[i] = mod62_mul_ypinv(wtab[i-1], w, wpinv, p);
	      winvtab[i] = mod62_mul_ypinv(winvtab[i-1], winv, winvpinv, p);
	    }
	  for (size_t i = 0; i < N / 2; i++)
	    {
	      wpinvtab[i] = mod62_ypinv(wtab[i], p, pinv);
	      winvpinvtab[i] = mod62_ypinv(winvtab[i], p, pinv);
	    }
	}
      else
	{
	  mod->wtab[lgN] = NULL;
	  mod->wpinvtab[lgN] = NULL;
	  mod->winvtab[lgN] = NULL;
	  mod->winvpinvtab[lgN] = NULL;
	}
    }
}


void fft62_mod_clear(fft62_mod_t* mod)
{
  for (unsigned lgN = 0; lgN <= FFT62_MAX_LGN; lgN++)
    if (mod->wtab[lgN] != NULL)
      {
	size_t N = (uint64_t) 1 << lgN;
	zz_free(mod->wtab[lgN], N/2 * sizeof(uint64_t));
	zz_free(mod->wpinvtab[lgN], N/2 * sizeof(uint64_t));
	zz_free(mod->winvtab[lgN], N/2 * sizeof(uint64_t));
	zz_free(mod->winvpinvtab[lgN], N/2 * sizeof(uint64_t));
      }
}



/*
  Given transform length N = 2^lgN, selects an array decomposition
  N = K * M, where
    K = number of rows = 2^lgK
    M = number of cols = 2^lgM
*/
void fft62_fft_array_params(unsigned* lgK, unsigned* lgM, unsigned lgN)
{
  *lgK = lgN / 2;
  *lgM = lgN - *lgK;
}


size_t fft62_next_size(size_t n, unsigned lgN)
{
  if (n == 0)
    n = 1;

  if (lgN <= FFT62_ARRAY_THRESHOLD)
    return n;

  unsigned lgK, lgM;
  fft62_fft_array_params(&lgK, &lgM, lgN);
  size_t k = ((n - 1) >> lgM) + 1;
  return fft62_next_size(k, lgK) << lgM;
}


static void
fft62_fft_layer_not_inplace(uint64_t* yp, uint64_t* xp, size_t blocks,
			    size_t size, uint64_t* wtab, uint64_t* wpinvtab,
			    uint64_t p)
{
  size /= 2;

  do
    {
      for (size_t j = 0; j < size; j++)
	{
	  uint64_t x0 = xp[j];              // x0 in [0, 2p)
	  uint64_t x1 = xp[j + size];       // x1 in [0, 2p)
	  yp[j] = mod63_add(x0, x1, 2*p);   // x0 + x1 in [0, 2p)
	  uint64_t t = x0 - x1 + 2*p;       // x0 - x1 in [0, 4p)
	  yp[j + size] = mod62_mul_ypinv_lazy(t, wtab[j], wpinvtab[j], p);
                                            // w*(x0 - x1) in [0, 2p)
	}

      xp += 2 * size;
      yp += 2 * size;
    }
  while (--blocks != 0);
}


// requires size divisible by 4
static void
fft62_fft_layer(uint64_t* xp, size_t blocks, size_t size,
		uint64_t* wtab, uint64_t* wpinvtab, uint64_t p)
{
  assert(size % 4 == 0);

  size /= 2;
  uint64_t* yp = xp + size;

  do
    {
      // 2-way unroll
      size_t j = size - 2;
      do
	{
	  uint64_t x0 = xp[j+1];              // x0 in [0, 2p)
	  uint64_t x1 = yp[j+1];              // x1 in [0, 2p)
	  uint64_t t0 = x0 - x1 + 2*p;        // x0 - x1 in [0, 4p)
	  yp[j+1] = mod62_mul_ypinv_lazy(t0, wtab[j+1], wpinvtab[j+1], p);
                                              // w*(x0 - x1) in [0, 2p)
	  xp[j+1] = mod63_add(x0, x1, 2*p);   // x0 + x1 in [0, 2p)

	  uint64_t x2 = xp[j];
	  uint64_t x3 = yp[j];
	  uint64_t t1 = x2 - x3 + 2*p;
	  yp[j] = mod62_mul_ypinv_lazy(t1, wtab[j], wpinvtab[j], p);
	  xp[j] = mod63_add(x2, x3, 2*p);
	}
      while (j -= 2);

      // last two butterflies

      uint64_t x0 = xp[1];
      uint64_t x1 = yp[1];
      xp[1] = mod63_add(x0, x1, 2*p);
      uint64_t t0 = x0 - x1 + 2*p;
      yp[1] = mod62_mul_ypinv_lazy(t0, wtab[1], wpinvtab[1], p);

      uint64_t x2 = xp[0];
      uint64_t x3 = yp[0];
      xp[0] = mod63_add(x2, x3, 2*p);
      yp[0] = mod63_sub(x2, x3, 2*p);

      xp += 2 * size;
      yp += 2 * size;
    }
  while (--blocks != 0);
}



static void
fft62_fft_last_two_layers(uint64_t* xp, size_t blocks,
			  uint64_t* wtab, uint64_t* wpinvtab, uint64_t p)
{
  // 4th root of unity
  uint64_t w = wtab[1];
  uint64_t wpinv = wpinvtab[1];

  do
    {
      uint64_t u0 = xp[0];
      uint64_t u1 = xp[1];
      uint64_t u2 = xp[2];
      uint64_t u3 = xp[3];

      uint64_t v0 = mod63_add(u0, u2, 2*p);
      uint64_t v2 = mod63_sub(u0, u2, 2*p);
      uint64_t v1 = mod63_add(u1, u3, 2*p);
      uint64_t v3 = mod62_mul_ypinv_lazy(u1 - u3 + 2*p, w, wpinv, p);

      xp[0] = mod63_add(v0, v1, 2*p);
      xp[1] = mod63_sub(v0, v1, 2*p);
      xp[2] = mod63_add(v2, v3, 2*p);
      xp[3] = mod63_sub(v2, v3, 2*p);

      xp += 4;
    }
  while (--blocks != 0);
}


void fft62_fft_base(uint64_t* yp, uint64_t* xp, unsigned lgN, fft62_mod_t* mod)
{
  if (lgN == 0)
    {
      yp[0] = xp[0];
      return;
    }

  uint64_t p = mod->p;

  if (lgN == 1)
    {
      uint64_t x0 = xp[0];
      uint64_t x1 = xp[1];
      yp[0] = mod63_add(x0, x1, 2*p);
      yp[1] = mod63_sub(x0, x1, 2*p);
      return;
    }

  uint64_t** wtab = mod->wtab;
  uint64_t** wpinvtab = mod->wpinvtab;

  if (lgN == 2)
    {
      if (yp != xp)
	yp[0] = xp[0], yp[1] = xp[1], yp[2] = xp[2], yp[3] = xp[3];
      fft62_fft_last_two_layers(yp, 1, wtab[2], wpinvtab[2], p);
      return;
    }

  unsigned j = lgN;
  size_t size = (size_t) 1 << lgN;
  size_t blocks = 1;

  if (yp != xp)
    {
      // do one layer out-of-place, to make next layer inplace
      fft62_fft_layer_not_inplace(yp, xp, blocks, size,
				  wtab[j], wpinvtab[j], p);
      j--, blocks <<= 1, size >>= 1;
    }

  for (; j > 2; j--, blocks <<= 1, size >>= 1)
    fft62_fft_layer(yp, blocks, size, wtab[j], wpinvtab[j], p);

  fft62_fft_last_two_layers(yp, blocks, wtab[2], wpinvtab[2], p);
}



void fft62_fft_short(uint64_t* yp, size_t yn, uint64_t* xp, size_t xn,
		     uint64_t* tp, unsigned lgN, fft62_mod_t* mod)
{
  size_t N = (size_t) 1 << lgN;

  assert(lgN <= FFT62_ARRAY_THRESHOLD);
  assert(xn >= 1 && xn <= N);
  assert(yn >= 1 && yn <= N);

  if (yn == N)
    {
      if (xn == N && lgN <= FFT62_CACHE_THRESHOLD)
	{
	  // no truncation
	  fft62_fft_base(yp, xp, lgN, mod);
	  return;
	}

      // use output as scratch space (so that the transform becomes inplace at
      // the next layer)
      tp = yp;
    }

  // divide-and-conquer algorithm

  size_t half = N >> 1;
  uint64_t p = mod->p;

  if (yn <= half)
    {
      if (xn <= half)
	{
	  fft62_fft_short(yp, yn, xp, xn, tp, lgN - 1, mod);
	}
      else
	{
	  xn -= half;

	  // (X, Y) -> X + Y
	  for (size_t j = 0; j < xn; j++)
	    tp[j] = mod63_add(xp[j], xp[j + half], 2*p);
	  if (tp != xp)
	    for (size_t j = xn; j < half; j++)
	      tp[j] = xp[j];

	  fft62_fft_short(yp, yn, tp, half, tp, lgN - 1, mod);
	}
    }
  else
    {
      yn -= half;
      tp += half;

      uint64_t* wtab = mod->wtab[lgN];
      uint64_t* wpinvtab = mod->wpinvtab[lgN];

      if (xn <= half)
	{
	  // X -> (X, w*X)
	  for (size_t j = 0; j < xn; j++)
	    tp[j] = mod62_mul_ypinv_lazy(xp[j], wtab[j], wpinvtab[j], p);

	  fft62_fft_short(yp, half, xp, xn, yp, lgN - 1, mod);
	  fft62_fft_short(yp + half, yn, tp, xn, tp, lgN - 1, mod);
	}
      else
	{
	  xn -= half;

	  // (X, Y) -> (X + Y, w*(X - Y))
	  for (size_t j = 0; j < xn; j++)
	    {
	      uint64_t x0 = xp[j];
	      uint64_t x1 = xp[j + half];
	      yp[j] = mod63_add(x0, x1, 2*p);
	      tp[j] = mod62_mul_ypinv_lazy(x0 - x1 + 2*p,
					   wtab[j], wpinvtab[j], p);
	    }
	  // X -> (X, w*X)
	  if (yp != xp)
	    for (size_t j = xn; j < half; j++)
	      yp[j] = xp[j];
	  for (size_t j = xn; j < half; j++)
	    tp[j] = mod62_mul_ypinv_lazy(xp[j], wtab[j], wpinvtab[j], p);

	  fft62_fft_short(yp, half, yp, half, yp, lgN - 1, mod);
	  fft62_fft_short(yp + half, yn, tp, half, tp, lgN - 1, mod);
	}
    }
}


// number of columns to grab at a time in array algorithms
#define GROUP 16

// amount of space to leave between adjacent columns
// (to reduce address conflicts in cache hits)
#define GAP 3


// helper functions for retrieving/storing columns out of an array

static void pull_columns(uint64_t* cp, uint64_t* xp, size_t cols, size_t rows,
			 size_t K, size_t M)
{
  for (size_t j = 0; j < rows; j++, cp++, xp += M)
    for (size_t i = 0; i < cols; i++)
      cp[i*K] = xp[i];
}


static void push_columns(uint64_t* xp, uint64_t* cp, size_t cols, size_t rows,
			 size_t K, size_t M)
{
  for (size_t j = 0; j < rows; j++, cp++, xp += M)
    for (size_t i = 0; i < cols; i++)
      xp[i] = cp[i*K];
}


// helper for fft62_fft_twisted()
static void
fft_twist(uint64_t* yp, uint64_t* xp, size_t xn, uint64_t z, unsigned lgH,
	  fft62_mod_t* mod)
{
  assert(z != 1);

  uint64_t p = mod->p;
  yp[0] = xp[0];

  if (z == 0 && lgH <= FFT62_ARRAY_THRESHOLD)
    {
      // get roots from table
      uint64_t* wtab = mod->wtab[lgH];
      uint64_t* wpinvtab = mod->wpinvtab[lgH];

      for (size_t i = 1; i < xn; i++)
     	yp[i] = mod62_mul_ypinv_lazy(xp[i], wtab[i], wpinvtab[i], p);
    }
  else
    {
      // compute roots on the fly
      if (z == 0)
	z = mod->rtab[lgH];

      uint64_t zpinv = mod62_ypinv(z, p, mod->pinv);
      uint64_t pinvb = mod->pinvb;
      // zi = z^i * B mod p, in [0, 2p)
      uint64_t zi = mod62_mul_ypinv_lazy(mod->B, z, zpinv, p);

      for (size_t i = 1; i < xn; i++)
	{
	  yp[i] = mod62_mul_pinvb_lazy(zi, xp[i], p, pinvb);
	  zi = mod62_mul_ypinv_lazy(zi, z, zpinv, p);
	}
    }
}



void fft62_fft_twisted(uint64_t* yp, size_t yn, uint64_t* xp, size_t xn,
		       unsigned lgN, uint64_t z, unsigned lgH,
		       fft62_mod_t* mod, int threads)
{
  size_t N = (size_t) 1 << lgN;

  if (lgN <= FFT62_ARRAY_THRESHOLD)
    {
      // small transform, use fft_short()
      uint64_t* tp = (yn == N) ? yp : zz_malloc(N * sizeof(uint64_t));

      // twist input coefficients if necessary
      if (z != 1)
	{
	  fft_twist(tp, xp, xn, z, lgH, mod);
	  xp = tp;
	}

      fft62_fft_short(yp, yn, xp, xn, tp, lgN, mod);
      if (tp != yp)
	zz_free(tp, N*sizeof(uint64_t));
      return;
    }

  // select array decomposition
  // K = number of rows = 2^lgK, M = number of columns = 2^lgM
  unsigned lgK, lgM;
  fft62_fft_array_params(&lgK, &lgM, lgN);
  size_t K = (size_t) 1 << lgK;
  size_t M = (size_t) 1 << lgM;

  // number of input and output rows
  assert(xn % M == 0);
  assert(yn % M == 0);
  size_t xr = xn >> lgM;
  size_t yr = yn >> lgM;

  // twist for column transforms = Z^M
  uint64_t zM = z;
  if (z > 1)
    {
      uint64_t p = mod->p;
      uint64_t pinv = mod->pinv;
      for (unsigned i = 0; i < lgM; i++)
	zM = mod62_mul_pinv_lazy2(zM, zM, p, pinv);
      zM = mod63_reduce2(zM, p);
    }
  unsigned lgHM = lgH - lgM;

  // transform columns
//#pragma omp parallel num_threads(threads)
  {
    uint64_t* cp = zz_malloc(GROUP * (K + GAP) * sizeof(uint64_t));

//#pragma omp for schedule(static)
    for (ptrdiff_t i = 0; i < M; i += GROUP)
      {
	// number of columns to handle on this iteration
	size_t cols = MIN(M - i, GROUP);

	// transfer columns to cp
	pull_columns(cp, xp + i, cols, xr, K + GAP, M);

	// transform columns
	uint64_t* tp = cp;
	for (size_t ii = 0; ii < cols; ii++, tp += K + GAP)
	  {
	    if (lgK <= FFT62_ARRAY_THRESHOLD)
	      {
		if (zM != 1)
		  fft_twist(tp, tp, xr, zM, lgHM, mod);
		fft62_fft_short(tp, yr, tp, xr, tp, lgK, mod);
	      }
	    else
	      fft62_fft_twisted(tp, yr, tp, xr, lgK, zM, lgHM, mod, 1);
	  }

	// transfer columns to yp
	push_columns(yp + i, cp, cols, yr, K + GAP, M);
      }

    zz_free(cp, GROUP * (K+GAP) * sizeof(uint64_t));
  }

  // distribute rows among available threads
//#pragma omp parallel for num_threads(threads) schedule(static)
  for (int t = 0; t < threads; t++)
    {
      size_t j_start = t * yr / threads;
      size_t j_end = (t + 1) * yr / threads;

      uint64_t* rp = yp + (j_start << lgM);

      // let w = basic root of order N.

      uint64_t p = mod->p;
      uint64_t pinv = mod->pinv;
      uint64_t* rtab = mod->rtab + lgN - lgK + 1;
      uint64_t* rinvtab = mod->rinvtab + lgN - lgK + 1;

      // Zwj = Z * w^(bit reversal of j), in [0, p)
      uint64_t Zwj = (z != 0) ? z : mod->rtab[lgH];
      for (unsigned k = 0; k < lgK; k++)
	if ((j_start >> k) & 1)
	  Zwj = mod62_mul_pinv_lazy4(Zwj, rtab[k], p, pinv);
      Zwj = mod62_reduce4(Zwj, p);

      for (size_t j = j_start; j < j_end; j++, rp += M)
	{
	  // transform row
	  fft62_fft_twisted(rp, M, rp, M, lgM, Zwj, 0, mod, 1);

	  // update Zwj
	  for (unsigned k = 0; k < lgK; k++)
	    {
	      if ((j >> k) & 1)
		Zwj = mod62_mul_pinv_lazy4(Zwj, rinvtab[k], p, pinv);
	      else
		{
		  Zwj = mod62_mul_pinv(Zwj, rtab[k], p, pinv);
		  break;
		}
	    }
	}
    }
}



// requires size divisible by 4
static void
fft62_ifft_layer(uint64_t* xp, size_t blocks, size_t size,
		 uint64_t* wtab, uint64_t* wpinvtab, uint64_t p)
{
  assert(size % 4 == 0);

  size /= 2;
  uint64_t* yp = xp + size;

  do
    {
      // 2-way unroll
      size_t j = size - 2;
      do
	{
	  uint64_t x1 = yp[j+1];            // x1 in [0, 4p)
	  uint64_t x3 = yp[j];
	  uint64_t t0 = mod62_mul_ypinv_lazy(x1, wtab[j+1], wpinvtab[j+1], p);
                                            // w*x1 in [0, 2p)
	  uint64_t t1 = mod62_mul_ypinv_lazy(x3, wtab[j], wpinvtab[j], p);

	  uint64_t x0 = mod63_reduce2(xp[j+1], 2*p);    // x0 in [0, 2p)
	  xp[j+1] = x0 + t0;                            // x0 + w*x1 in [0, 4p)
	  yp[j+1] = x0 - t0 + 2*p;                      // x0 - w*x1 in [0, 4p)

	  uint64_t x2 = mod63_reduce2(xp[j], 2*p);
	  xp[j] = x2 + t1;
	  yp[j] = x2 - t1 + 2*p;
	}
      while (j -= 2);

      // last two butterflies

      uint64_t x0 = xp[1];
      uint64_t x1 = yp[1];
      uint64_t t0 = mod62_mul_ypinv_lazy(x1, wtab[1], wpinvtab[1], p);
      x0 = mod63_reduce2(x0, 2*p);
      xp[1] = x0 + t0;
      yp[1] = x0 - t0 + 2*p;

      uint64_t x2 = mod63_reduce2(xp[0], 2*p);
      uint64_t x3 = mod63_reduce2(yp[0], 2*p);
      xp[0] = x2 + x3;
      yp[0] = x2 - x3 + 2*p;

      xp += 2 * size;
      yp += 2 * size;
    }
  while (--blocks != 0);
}


static void
fft62_ifft_first_two_layers(uint64_t* yp, uint64_t* xp, size_t blocks,
			    uint64_t* wtab, uint64_t* wpinvtab, uint64_t p)
{
  // 4th root of unity
  uint64_t w = wtab[1];
  uint64_t wpinv = wpinvtab[1];

  do
    {
      uint64_t u0 = mod63_reduce2(xp[0], 2*p);
      uint64_t u1 = mod63_reduce2(xp[1], 2*p);
      uint64_t u2 = mod63_reduce2(xp[2], 2*p);
      uint64_t u3 = mod63_reduce2(xp[3], 2*p);

      uint64_t v0 = mod63_add(u0, u1, 2*p);
      uint64_t v1 = mod63_sub(u0, u1, 2*p);
      uint64_t v2 = mod63_add(u2, u3, 2*p);
      uint64_t v3 = mod62_mul_ypinv_lazy(u2 - u3 + 2*p, w, wpinv, p);

      yp[0] = v0 + v2;
      yp[2] = v0 - v2 + 2*p;
      yp[1] = v1 + v3;
      yp[3] = v1 - v3 + 2*p;

      xp += 4;
      yp += 4;
    }
  while (--blocks != 0);
}



void fft62_ifft_base(uint64_t* yp, uint64_t* xp, unsigned lgN, fft62_mod_t* mod)
{
  if (lgN == 0)
    {
      yp[0] = xp[0];
      return;
    }

  uint64_t p = mod->p;

  if (lgN == 1)
    {
      uint64_t x0 = mod63_reduce2(xp[0], 2*p);
      uint64_t x1 = mod63_reduce2(xp[1], 2*p);
      yp[0] = x0 + x1;
      yp[1] = x0 - x1 + 2*p;
      return;
    }

  uint64_t** wtab = mod->winvtab;
  uint64_t** wpinvtab = mod->winvpinvtab;

  size_t blocks = (size_t) 1 << (lgN - 2);
  fft62_ifft_first_two_layers(yp, xp, blocks, wtab[2], wpinvtab[2], p);
  blocks >>= 1;

  size_t size = 8;
  for (unsigned j = 3; j <= lgN; j++, blocks >>= 1, size <<= 1)
    fft62_ifft_layer(yp, blocks, size, wtab[j], wpinvtab[j], p);
}



void fft62_ifft_short1(uint64_t* yp, size_t yn, uint64_t* xp,
		       uint64_t* tp, unsigned lgN, fft62_mod_t* mod)
{
  size_t N = (size_t) 1 << lgN;

  assert(lgN <= FFT62_ARRAY_THRESHOLD);
  assert(1 <= yn && yn <= N);

  if (yn == N && lgN <= FFT62_CACHE_THRESHOLD)
    {
      // no truncation
      fft62_ifft_base(yp, xp, lgN, mod);
      return;
    }

  // divide-and-conquer algorithm

  size_t half = N >> 1;
  uint64_t p = mod->p;

  if (yn <= half)
    {
      // X -> 2X
      for (size_t j = 0; j < yn; j++)
      	tp[j] = 2 * mod63_reduce2(xp[j], 2*p);

      fft62_ifft_short1(yp, yn, tp, tp, lgN - 1, mod);
    }
  else
    {
      fft62_ifft_short1(yp, half, xp, yp, lgN - 1, mod);

      yn -= half;
      tp += half;

      uint64_t* wtab = mod->wtab[lgN];
      uint64_t* wpinvtab = mod->wpinvtab[lgN];

      if (tp != xp + half)
	for (size_t j = 0; j < yn; j++)
      	  tp[j] = xp[j + half];
      // X -> (2X, w*X)
      for (size_t j = yn; j < half; j++)
	{
	  uint64_t x0 = yp[j];
	  tp[j] = mod62_mul_ypinv_lazy(x0, wtab[j], wpinvtab[j], p);
	  yp[j] = 2 * mod63_reduce2(x0, 2*p);
	}

      fft62_ifft_short2(tp, yn, tp, tp, lgN - 1, mod);

      wtab = mod->winvtab[lgN];
      wpinvtab = mod->winvpinvtab[lgN];

      // (X, Y) -> (X + Y/w, X - Y/w)
      for (size_t j = 0; j < yn; j++)
	{
	  uint64_t x0 = mod63_reduce2(yp[j], 2*p);
	  uint64_t x1 = mod62_mul_ypinv_lazy(tp[j], wtab[j], wpinvtab[j], p);
	  yp[j] = x0 + x1;
	  yp[j + half] = x0 - x1 + 2*p;
	}
    }
}



void fft62_ifft_short2(uint64_t* yp, size_t yn, uint64_t* xp,
		       uint64_t* tp, unsigned lgN, fft62_mod_t* mod)
{
  size_t N = (size_t) 1 << lgN;

  assert(lgN <= FFT62_ARRAY_THRESHOLD);
  assert(1 <= yn && yn <= N);

  if (yn == N && lgN <= FFT62_CACHE_THRESHOLD)
    {
      // no truncation
      fft62_ifft_base(yp, xp, lgN, mod);
      return;
    }

  // divide-and-conquer algorithm

  size_t half = N >> 1;
  uint64_t p = mod->p;

  if (yn <= half)
    {
      // X -> 2X
      for (size_t j = 0; j < yn; j++)
     	tp[j] = 2 * mod63_reduce2(xp[j], 2*p);
      // (X, Y) -> X + Y
      for (size_t j = yn; j < half; j++)
	tp[j] = mod64_add(xp[j], xp[j + half], 4*p);

      fft62_ifft_short2(yp, yn, tp, tp, lgN - 1, mod);

      // (X, Y) -> X - Y
      for (size_t j = 0; j < yn; j++)
	yp[j] = mod64_sub(yp[j], xp[j + half], 4*p);
    }
  else
    {
      fft62_ifft_short1(yp, half, xp, yp, lgN - 1, mod);

      yn -= half;
      tp += half;

      uint64_t* wtab = mod->wtab[lgN];
      uint64_t* wpinvtab = mod->wpinvtab[lgN];

      if (tp != xp + half)
	for (size_t j = 0; j < yn; j++)
      	  tp[j] = xp[j + half];
      // (X, Y) -> (2X - Y, w*(X - Y))
      for (size_t j = yn; j < half; j++)
	{
	  uint64_t x0 = yp[j];
	  uint64_t x1 = xp[j + half];
	  uint64_t u = mod64_sub(x0, x1, 4*p);
	  tp[j] = mod62_mul_ypinv_lazy(u, wtab[j], wpinvtab[j], p);
	  yp[j] = mod64_add(x0, u, 4*p);
	}

      fft62_ifft_short2(tp, yn, tp, tp, lgN - 1, mod);

      wtab = mod->winvtab[lgN];
      wpinvtab = mod->winvpinvtab[lgN];

      // (X, Y) -> (X + Y/w, X - Y/w)
      for (size_t j = 0; j < yn; j++)
	{
	  uint64_t x0 = mod63_reduce2(yp[j], 2*p);
	  uint64_t x1 = mod62_mul_ypinv_lazy(tp[j], wtab[j], wpinvtab[j], p);
	  yp[j] = x0 + x1;
	  yp[j + half] = x0 - x1 + 2*p;
	}
    }
}



// helper for fft62_ifft_twisted()
// (similar to fft_twist(), but always operates inplace, and inputs and outputs
// are in [0, 4p) instead of [0, 2p))
static void
ifft_twist(uint64_t* xp, size_t xn, uint64_t z, unsigned lgH,
	   fft62_mod_t* mod)
{
  assert(z != 1);

  uint64_t p = mod->p;

  if (z == 0 && lgH <= FFT62_ARRAY_THRESHOLD)
    {
      // get roots from table
      uint64_t* wtab = mod->winvtab[lgH];
      uint64_t* wpinvtab = mod->winvpinvtab[lgH];

      for (size_t i = 1; i < xn; i++)
	xp[i] = mod62_mul_ypinv_lazy(xp[i], wtab[i], wpinvtab[i], p);
    }
  else
    {
      // compute roots on the fly
      if (z == 0)
	z = mod->rinvtab[lgH];

      uint64_t zpinv = mod62_ypinv(z, p, mod->pinv);
      uint64_t pinvb = mod->pinvb;
      // zi = z^i * B mod p, in [0, 2p)
      uint64_t zi = mod62_mul_ypinv_lazy(mod->B, z, zpinv, p);

      for (size_t i = 1; i < xn; i++)
	{
	  xp[i] = mod62_mul_pinvb_lazy(zi, xp[i], p, pinvb);
	  zi = mod62_mul_ypinv_lazy(zi, z, zpinv, p);
	}
    }
}



void fft62_ifft_twisted(uint64_t* yp, size_t yn, uint64_t* xp,
			unsigned lgN, uint64_t z, unsigned lgH,
			fft62_mod_t* mod, int threads)
{
  size_t N = (size_t) 1 << lgN;

  if (lgN <= FFT62_ARRAY_THRESHOLD)
    {
      // small transform, use ifft_short()
      uint64_t* tp = (yn == N) ? yp : zz_malloc(N * sizeof(uint64_t));
      fft62_ifft_short1(yp, yn, xp, tp, lgN, mod);
      if (tp != yp)
	zz_free(tp, N * sizeof(uint64_t));

      // untwist output coefficients if necessary
      if (z != 1)
	ifft_twist(yp, yn, z, lgH, mod);

      return;
    }

  // select array decomposition
  // K = number of rows = 2^lgK, M = number of columns = 2^lgM
  unsigned lgK, lgM;
  fft62_fft_array_params(&lgK, &lgM, lgN);
  size_t K = (size_t) 1 << lgK;
  size_t M = (size_t) 1 << lgM;

  // number of input and output rows
  assert(yn % M == 0);
  size_t yr = yn >> lgM;

  // distribute rows among available threads
//#pragma omp parallel for num_threads(threads) schedule(static)
  for (int t = 0; t < threads; t++)
    {
      size_t j_start = t * yr / threads;
      size_t j_end = (t + 1) * yr / threads;

      uint64_t* rp = xp + (j_start << lgM);
      uint64_t* sp = yp + (j_start << lgM);

      // let w = basic root of order N.

      uint64_t p = mod->p;
      uint64_t pinv = mod->pinv;
      uint64_t* rtab = mod->rtab + lgN - lgK + 1;
      uint64_t* rinvtab = mod->rinvtab + lgN - lgK + 1;

      // Zwj = (Z * w^(bit reversal of j))^(-1), in [0, p)
      uint64_t Zwj = (z != 0) ? z : mod->rinvtab[lgH];
      for (unsigned k = 0; k < lgK; k++)
	if ((j_start >> k) & 1)
	  Zwj = mod62_mul_pinv_lazy4(Zwj, rinvtab[k], p, pinv);
      Zwj = mod62_reduce4(Zwj, p);

      for (size_t j = j_start; j < j_end; j++, rp += M, sp += M)
	{
	  // transform row
	  fft62_ifft_twisted(sp, M, rp, lgM, Zwj, 0, mod, 1);

	  // update Zwj
	  for (unsigned k = 0; k < lgK; k++)
	    {
	      if ((j >> k) & 1)
		Zwj = mod62_mul_pinv_lazy4(Zwj, rtab[k], p, pinv);
	      else
		{
		  Zwj = mod62_mul_pinv(Zwj, rinvtab[k], p, pinv);
		  break;
		}
	    }
	}
    }

  // twist for column transforms = Z^M
  uint64_t zM = z;
  if (z > 1)
    {
      uint64_t p = mod->p;
      uint64_t pinv = mod->pinv;
      for (unsigned i = 0; i < lgM; i++)
	zM = mod62_mul_pinv_lazy2(zM, zM, p, pinv);
      zM = mod63_reduce2(zM, p);
    }
  unsigned lgHM = lgH - lgM;

  // transform columns
//#pragma omp parallel num_threads(threads)
  {
    uint64_t* cp = zz_malloc(GROUP * (K + GAP) * sizeof(uint64_t));

//#pragma omp for schedule(static)
    for (ptrdiff_t i = 0; i < M; i += GROUP)
      {
	// number of columns to handle on this iteration
	size_t cols = MIN(M - i, GROUP);

	// transfer columns to cp
	pull_columns(cp, yp + i, cols, yr, K + GAP, M);

	// transform columns
	uint64_t* tp = cp;
	for (size_t ii = 0; ii < cols; ii++, tp += K + GAP)
	  {
	    if (lgK <= FFT62_ARRAY_THRESHOLD)
	      {
		fft62_ifft_short1(tp, yr, tp, tp, lgK, mod);
		if (zM != 1)
		  ifft_twist(tp, yr, zM, lgHM, mod);
	      }
	    else
	      fft62_ifft_twisted(tp, yr, tp, lgK, zM, lgHM, mod, 1);
	  }

	// transfer columns to yp
	push_columns(yp + i, cp, cols, yr, K + GAP, M);
      }

    zz_free(cp, GROUP * (K + GAP) * sizeof(uint64_t));
  }
}
