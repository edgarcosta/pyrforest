/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include <assert.h>
#include <stdint.h>
#include <stddef.h>
#include "mpzfft_moduli.h"
#include "zzcrt.h"
#include "recompose.h"
#include "zzmisc.h"
#include "zzmem.h"

void zz_crt_recompose(mp_limb_t* rp, size_t rn, unsigned r,
		      uint64_t** up, size_t un,
		      zz_moduli_t* moduli, unsigned num_primes, int threads)
{
  // chunk size, want everything to safely fit into L1
  // (must be divisible by 64)
  const size_t size = 512;

  size_t chunks = (un + size - 1) / size;

  // size of output region corresponding to each chunk
  size_t len = size * r / 64;

  // t = number of coefficients to process directly into output for each chunk
  size_t t;
  if (size * r < 62 * num_primes + 1)
    t = 0;
  else
    {
      t = size - (62 * num_primes / r);
      assert(t >= 1);
    }

  // sp_size = number of limbs needed for recomposition of last few coefficients
  // for each chunk
  size_t sp_size = ((size - 1) * r + 62 * num_primes) / 64 - (t * r) / 64 + 1;
  assert((ptrdiff_t) sp_size >= 0);
  if (sp_size == 0)
    sp_size = 1;
  mp_limb_t* sp = zz_malloc(chunks * sp_size * sizeof(mp_limb_t));

//#pragma omp parallel num_threads(threads)
  {
    mp_limb_t* temp = zz_malloc(num_primes * size * sizeof(mp_limb_t));
    uint64_t* up_thread[num_primes];

  // process "up" in chunks
//#pragma omp for schedule(static)
    for (ptrdiff_t chunk = 0; chunk < chunks; chunk++)
      {
	// this iteration handles inputs i0 <= i < i1
	size_t i0 = chunk * size;
	size_t i1 = MIN(i0 + size, un);

	// write output starting at offset out0
	size_t out0 = i0 * r / 64;
	if (out0 >= rn)
	  continue;     // this chunk does not contribute to output
	size_t out1 = MIN(rn, out0 + len);

	// CRT for this chunk, write coefficients to temp
	for (unsigned j = 0; j < num_primes; j++)
	  up_thread[j] = up[j] + i0;
	zz_crt(temp, up_thread, i1 - i0, moduli, num_primes);

	// recompose coefficients i0 <= i < i0 + t directly to output
	if (t > 0)
	  zz_recompose(rp + out0, out1 - out0, temp, MIN(t, i1 - i0),
		       num_primes, r, 0);
	else
	  mpn_zero(rp + out0, out1 - out0);

	// recompose coefficients i0 + t <= i < i1 to sp, to be added in later
	mp_limb_t* sp_dest = sp + chunk * sp_size;
	if (t < i1 - i0)
	  zz_recompose(sp_dest, sp_size, temp + t * num_primes,
		       i1 - i0 - t, num_primes, r, (t * r) % 64);
	else
	  mpn_zero(sp_dest, sp_size);

	// incorporate sign of main part into sp part
	if (out1 < rn)
	  {
	    size_t k = len - t * r / 64;
	    mpn_sub_1(sp_dest + k, sp_dest + k, sp_size - k,
		      rp[out1 - 1] >> 63);
	  }
      }

    zz_free(temp, num_primes * size *sizeof(mp_limb_t));
  }

  // zero-pad up to length rn
  size_t written = chunks * len;
  if (written < rn)
    mpn_zero(rp + written, rn - written);
  
  // add in missing pieces from sp
  for (size_t chunk = 0; chunk < chunks; chunk++)
    {
      mp_limb_t* src = sp + chunk * sp_size;
      mp_limb_t* dest = rp + chunk * len + t * r / 64;

      if (dest >= rp + rn)
	break;

      if (sp_size >= rp + rn - dest)
	{
	  mpn_add_n(dest, dest, src, rp + rn - dest);
	  continue;
	}
      mp_limb_signed_t cy = mpn_add_n(dest, dest, src, sp_size);
      cy -= src[sp_size - 1] >> 63;

      if (chunk == chunks - 1)
	{
	  // propagate carry/borrow to end of output buffer
	  if (cy >= 0)
	    mpn_add_1(dest + sp_size, dest + sp_size,
		      (rp + rn) - (dest + sp_size), cy);
	  else
	    mpn_sub_1(dest + sp_size, dest + sp_size,
		      (rp + rn) - (dest + sp_size), -cy);
	  break;
	}

      // propagate carry/borrow into next piece of sp
      if (len > sp_size)
	{
	  size_t k = MIN(len - sp_size, (rp + rn) - (dest + sp_size));
	  if (cy >= 0)
	    {
	      cy = mpn_add_1(dest + sp_size, dest + sp_size, k, cy);
	      mpn_add_1(src + sp_size, src + sp_size, sp_size, cy);
	    }
	  else
	    {
	      cy = mpn_sub_1(dest + sp_size, dest + sp_size, k, -cy);
	      mpn_sub_1(src + sp_size, src + sp_size, sp_size, cy);
	    }
	}
      else
	{
	  if (cy >= 0)
	    mpn_add_1(src + 2*sp_size - len, src + 2*sp_size - len, len, cy);
	  else
	    mpn_sub_1(src + 2*sp_size - len, src + 2*sp_size - len, len, -cy);
	}
    }

  zz_free(sp, chunks * sp_size * sizeof(mp_limb_t));
}
