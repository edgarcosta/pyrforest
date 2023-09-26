/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ARITH128_H
#define ARITH128_H

#include <stdint.h>


/*
  Functions for 128-bit arithmetic (all parameters are uint64_t):

  MUL_WIDE(z1, z0, x, y):
    computes z0 and z1 so that 2^64*z1 + z0 = x * y

  MUL_HI(x, y):
    like MUL_WIDE, but returns only z1

  ADD_WIDE(z1, z0, x1, x0, y1, y0):
    computes z0 and z1 so that
        (2^64*x1 + x0) + (2^64*y1 + y0) = 2^64*z1 + z0  mod 2^128

  DIV_WIDE(q, r, x1, x0, y):
    computes q and r so that (2^64*x1 + x0) = q * y + r,
    assuming the quotient fits into 64 bits

  These are implemented either by:
    * inline assembly if Z_USE_INLINE_ASM is defined (default), or
    * native 128-bit types if Z_USE_128_BIT_TYPES is defined.

  note: the inline assembly version seems to generate better assembly for
  ADD_WIDE than gcc's native 128-bit capability, at least up to gcc 4.7.2.
  See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=51838
*/


// Warning: Z_USE_INLINE_ASM currently fails when used with clang
#if !defined(Z_USE_128_BIT_TYPES) && !defined(Z_USE_INLINE_ASM)
#define Z_USE_128_BIT_TYPES
#endif


#ifdef Z_USE_128_BIT_TYPES

typedef int int128_t __attribute__((mode(TI)));
typedef unsigned int uint128_t __attribute__((mode(TI)));

#define MUL_WIDE(z1, z0, x, y)			\
  do {						\
    uint128_t __z = (uint128_t) (x) * (y);	\
    (z0) = (uint64_t) __z;			\
    (z1) = __z >> 64;				\
  } while (0)

#define ADD_WIDE(z1, z0, x1, x0, y1, y0)		\
  do {							\
    uint128_t __x = ((uint128_t) x1 << 64) + x0;	\
    uint128_t __y = ((uint128_t) y1 << 64) + y0;	\
    uint128_t __z = __x + __y;				\
    (z0) = (uint64_t) __z;				\
    (z1) = __z >> 64;					\
  } while (0)

#define DIV_WIDE(q, r, x1, x0, y)			\
  do {							\
    uint128_t __x = ((uint128_t) (x1) << 64) + (x0);	\
    uint64_t __y = (y);					\
    (q) = __x / __y;					\
    (r) = __x % __y;					\
  } while (0)

#endif


#ifdef Z_USE_INLINE_ASM

#define MUL_WIDE(z1, z0, x, y) __asm__("mulq %3" : "=a" (z0), "=d" (z1) : "%0" ((uint64_t)(x)), "rm" ((uint64_t)(y)))
#define DIV_WIDE(q, r, x1, x0, y) __asm__("divq %4" : "=a" (q), "=d" (r) : "0" ((uint64_t)(x0)), "1" ((uint64_t)(x1)), "rm" ((uint64_t)(y)))
#define ADD_WIDE(z1, z0, x1, x0, y1, y0) __asm__ ("addq %5, %1\n\tadcq %3, %0" : "=r,m" (z1), "=&r,m" (z0) : "0,0" ((uint64_t)(x1)), "rme,re" ((uint64_t)(y1)), "1,1" ((uint64_t)(x0)), "rme,re" ((uint64_t)(y0)))

#endif


static inline uint64_t MUL_HI(uint64_t x, uint64_t y)
{
  uint64_t z1, __attribute__((unused)) z0;
  MUL_WIDE(z1, z0, x, y);
  return z1;
}


#endif
