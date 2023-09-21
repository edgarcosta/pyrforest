/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef FFT62_FFT62_H
#define FFT62_FFT62_H

#include <stdlib.h>
#include <stdint.h>


/*
  FFT transform lengths are of the form N = 2^lgN with
  0 <= lgN <= FFT62_MAX_LGN.

  All FFT primes must satisfy p = 1 mod 2^FFT62_MAX_LGN.
*/

#define FFT62_MAX_LGN 44


// some primes p = 1 mod 2^44:
#define FFT62_PRIME1 ((uint64_t) 0x3fffc00000000001ULL)
#define FFT62_PRIME2 ((uint64_t) 0x3ffdf00000000001ULL)
#define FFT62_PRIME3 ((uint64_t) 0x3ffd900000000001ULL)
#define FFT62_PRIME4 ((uint64_t) 0x3ffd500000000001ULL)
#define FFT62_PRIME5 ((uint64_t) 0x3ffc100000000001ULL)
#define FFT62_PRIME6 ((uint64_t) 0x3ffbb00000000001ULL)
#define FFT62_PRIME7 ((uint64_t) 0x3ffaf00000000001ULL)
#define FFT62_PRIME8 ((uint64_t) 0x3ffac00000000001ULL)



/*
  The DFT is defined as follows.

  Let the input sequence be a_0, ..., a_{N-1}.

  Let w = standard primitive N-th root of 1, i.e. w = g^(2^FFT62_MAX_LGN / N),
  where g = some fixed element of Z/pZ of order 2^FFT62_MAX_LGN.

  Let Z = an element of (Z/pZ)^* (twisting parameter).

  Then the output sequence is
    b_j = \sum_{0 <= i < N} Z^i a_i w^(ij'), for 0 <= j < N,
  where j' is the length-lgN bit-reversal of j.

  Some of the FFT routines can operate on truncated sequences of certain
  "admissible" sizes. A size parameter n is admissible if 1 <= n <= N, and n is
  divisible by a certain power of 2. The precise power depends on the recursive
  array decomposition of the FFT. The smallest admissible n' >= n can be
  obtained via fft62_next_size().
*/



// transform sizes >= 2^FFT62_ARRAY_THRESHOLD use array decomposition algorithm
// transform sizes >= 2^FFT62_CACHE_THRESHOLD are split into two subtransforms

#define FFT62_ARRAY_THRESHOLD_DEFAULT 16
#define FFT62_CACHE_THRESHOLD_DEFAULT 13

#if TEST
// in the test suite, we want to be able to modify the thresholds to test
// certain algorithms

extern unsigned FFT62_ARRAY_THRESHOLD;
extern unsigned FFT62_CACHE_THRESHOLD;

#else
// in a regular build, want these to be compile-time constants

#define FFT62_ARRAY_THRESHOLD FFT62_ARRAY_THRESHOLD_DEFAULT
#define FFT62_CACHE_THRESHOLD FFT62_CACHE_THRESHOLD_DEFAULT

#endif



typedef struct
{
  uint64_t p, pinv, pinvb;

  uint64_t B;     // 2^64 mod p

  // g = primitive (2^FFT62_MAX_LGN)-th root mod p
  // ginv = 1/g mod p
  uint64_t g, ginv;

  // rtab[lgN] = basic root of order 2^lgN
  uint64_t rtab[FFT62_MAX_LGN + 1];
  // rinvtab[j] = 1 / rtab[j] mod p
  uint64_t rinvtab[FFT62_MAX_LGN + 1];

  /*
    wtab is an array of root tables. For a given N = 2^lgN, if
    1 <= lgN <= FFT62_ARRAY_THRESHOLD, then wtab[lgN] points to a table of
    length N/2, with the roots {1, w, w^2, ..., w^(N/2-1)}, where
    w = g^(2^FFT62_MAX_LGN / N). All roots are in [0, p).
    Otherwise wtab[lgN] == NULL.
  */
  uint64_t* wtab[FFT62_MAX_LGN + 1];

  // wpinvtab[j] = floor(2^64 * wtab[j] / p), for use with mod62_mul_ypinv()
  uint64_t* wpinvtab[FFT62_MAX_LGN + 1];

  // same as wtab and wpinvtab, but with w replaced by 1/w
  uint64_t* winvtab[FFT62_MAX_LGN + 1];
  uint64_t* winvpinvtab[FFT62_MAX_LGN + 1];
}
fft62_mod_t;


// initialises mod for use with p
void fft62_mod_init(fft62_mod_t* mod, uint64_t p);


// destroys mod
void fft62_mod_clear(fft62_mod_t* mod);


// ceil(log2(n)) if n >= 1, or 0 if n == 0
unsigned fft62_log2(size_t n);


// returns smallest admissible size n' >= n (see above)
size_t fft62_next_size(size_t n, unsigned lgN);


/*
  Truncated FFT interface is as follows:

  xn and yn must be admissible sizes for N.

  Input in xp[] is a_0, a_1, ..., a_{xn-1}. Assumes a_i = 0 for xn <= i < N.

  Output in yp[] is b_0, ..., b_{yn-1}, i.e. only first yn outputs are computed.

  Twisting parameter Z is described by z and lgH. If z == 0, then Z = basic
  2^lgH-th root of 1, and must have lgH >= lgN + 1. If z != 0, then Z = z
  (and lgH is ignored).

  The buffers {xp,xn} and {yp,yn} may overlap, but only if xp == yp.

  Inputs are in [0, 2p), outputs are in [0, 2p).

  threads = number of OpenMP threads to use.
*/

// implements full truncated FFT interface.
void fft62_fft_twisted(uint64_t* yp, size_t yn, uint64_t* xp, size_t xn,
		       unsigned lgN, uint64_t z, unsigned lgH,
		       fft62_mod_t* mod, int threads);

// special case of fft62_fft_twisted() for Z = 1
static inline
void fft62_fft(uint64_t* yp, size_t yn, uint64_t* xp, size_t xn, unsigned lgN,
	       fft62_mod_t* mod, int threads)
{
  fft62_fft_twisted(yp, yn, xp, xn, lgN, 1, 0, mod, threads);
}


/*
  Implements truncated FFT interface, with these restrictions:
  * Always Z = 1 and threads = 1.
  * Must have lgN <= FFT62_ARRAY_THRESHOLD.
  * Scratch space must also be provided in {tp,N}. The buffer {tp,N} may overlap
    with {xp,xn} (and/or {yp,yn}), but only if tp == xp (and/or tp == yp).
*/
void fft62_fft_short(uint64_t* yp, size_t yn, uint64_t* xp, size_t xn,
		     uint64_t* tp, unsigned lgN, fft62_mod_t* mod);


/*
  Implements truncated FFT interface, with these restrictions:
  * Always Z = 1 and threads = 1.
  * Must have lgN <= FFT62_ARRAY_THRESHOLD.
  * Always non-truncated, i.e. xn == yn == N.
*/
void fft62_fft_base(uint64_t* yp, uint64_t* xp, unsigned lgN, fft62_mod_t* mod);


/*
  Inverse truncated FFT interface is as follows.

  xn and yn must be admissible sizes for N, with yn <= xn.

  Input in xp[] is b_0, b_1, ..., b_{yn-1}, N*a_{yn}, ..., N*a_{xn-1}.

  Assumes a_i = 0 for xn <= i < N.

  Output in yp[] is N*a_0, ..., N*a_{yn-1}.

  Twisting parameter Z is described by z and lgH. If z == 0, then Z = basic
  2^lgH-th root of 1, and must have lgH >= lgN + 1. If z != 0, then Z = z^(-1)
  (and lgH is ignored).

  The buffers {xp,xn} and {yp,yn} may overlap, but only if xp == yp.

  Inputs are in [0, 4p), outputs are in [0, 4p).

  threads = number of OpenMP threads to use.

  (note: no function actually implements this interface in full generality!
  This is because it is tricky (and not that useful) to implement the twisting
  parameter when xn != yn.)
*/

/*
  Implements truncated inverse FFT interface, with the restriction that
  xn == yn.
*/
void fft62_ifft_twisted(uint64_t* yp, size_t yn, uint64_t* xp, unsigned lgN,
			uint64_t z, unsigned lgH,
			fft62_mod_t* mod, int threads);

// special case of fft62_ifft_twisted() for Z = 1
static inline
void fft62_ifft(uint64_t* yp, size_t yn, uint64_t* xp, unsigned lgN,
		fft62_mod_t* mod, int threads)
{
  fft62_ifft_twisted(yp, yn, xp, lgN, 1, 0, mod, threads);
}


/*
  Implements truncated inverse FFT interface, with these restrictions:
  * Always Z = 1 and threads = 1.
  * Always xn == yn for short1(), and always xn == N for short2().
  * Must have lgN <= FFT62_ARRAY_THRESHOLD.
  * Scratch space must also be provided in {tp,N}. The buffer {tp,N} may overlap
    with {xp,xn} (and/or {yp,yn}), but only if tp == xp (and/or tp == yp).
*/
void fft62_ifft_short1(uint64_t* yp, size_t yn, uint64_t* xp,
		       uint64_t* tp, unsigned lgN, fft62_mod_t* mod);
void fft62_ifft_short2(uint64_t* yp, size_t yn, uint64_t* xp,
		       uint64_t* tp, unsigned lgN, fft62_mod_t* mod);

/*
  Implements truncated inverse FFT interface, with these restrictions:
  * Always Z = 1 and threads = 1.
  * Always xn == yn == N.
  * Must have lgN <= FFT62_ARRAY_THRESHOLD.
*/
void fft62_ifft_base(uint64_t* yp, uint64_t* xp, unsigned lgN,
		     fft62_mod_t* mod);


#endif
