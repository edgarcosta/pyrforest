/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef FFT62_MOD62_H
#define FFT62_MOD62_H


#include <assert.h>
#include "arith128.h"


// in this file, B = 2^64


/*
  modular/addition subtraction for moduli up to 63 or 64 bits
*/

// m in (0, B/2), x, y in [0, m)
// returns x + y mod m, in [0, m)
static inline uint64_t mod63_add(uint64_t x, uint64_t y, uint64_t m)
{
  uint64_t z = x + y;
  z -= (z >= m) ? m : 0;
  return z;
}


// m in (0, B/2), x, y in [0, m)
// returns x - y mod m, in [0, m)
static inline uint64_t mod63_sub(uint64_t x, uint64_t y, uint64_t m)
{
  uint64_t z = x - y;
  z += ((int64_t) z < 0) ? m : 0;
  return z;
}


// m in (0, B), x, y in [0, m)
// returns x + y mod m, in [0, m)
static inline uint64_t mod64_add(uint64_t x, uint64_t y, uint64_t m)
{
  y = m - y;
  uint64_t z = x - y;
  z += (x < y) ? m : 0;
  return z;
}


// m in (0, B), x, y in [0, m)
// returns x - y mod m, in [0, m)
static inline uint64_t mod64_sub(uint64_t x, uint64_t y, uint64_t m)
{
  uint64_t z = x - y;
  z += (x < y) ? m : 0;
  return z;
}


// m in (0, B/2), x in [0, 2m)
// return x mod m, in [0, m)
static inline uint64_t mod63_reduce2(uint64_t x, uint64_t m)
{
  return x - ((x >= m) ? m : 0);
}


// m in (0, B/4), x in [0, 4m)
// return x mod m, in [0, m)
static inline uint64_t mod62_reduce4(uint64_t x, uint64_t m)
{
  return mod63_reduce2(mod63_reduce2(x, 2*m), m);
}



// *****************************************************************************
// 62-bit modular arithmetic

// in this section, p is odd and satisfies B/8 < p < B/4


// checks conditions on p
static inline int mod62_valid(uint64_t p)
{
  return (p & 1) && (p > 0x2000000000000000ULL) && (p < 0x4000000000000000ULL);
}


// x in [0, 2p)
// returns x / 2 mod p, in [0, 2p)
static inline uint64_t mod62_div_by2(uint64_t x, uint64_t p)
{
  return (x + ((x & 1) ? p : 0)) >> 1;
}


// x * y mod p, stupid hardware division
static inline uint64_t mod62_mul(uint64_t x, uint64_t y, uint64_t p)
{
  uint64_t z0, z1, __attribute__((unused)) q, r; // silence pedantic "unused-but-set-variable"
  MUL_WIDE(z1, z0, x, y);
  DIV_WIDE(q, r, z1, z0, p);
  return r;
}


// floor(B^2 / 4p) - B, in [0, B)
static inline uint64_t mod62_pinv(uint64_t p)
{
  assert(mod62_valid(p));

  uint64_t q, __attribute__((unused)) r; // silence pedantic "unused-but-set-variable"
  DIV_WIDE(q, r, -4*p, 0, 4*p);
  return q;
}



// 1 / p mod B, in [0, B)
static inline uint64_t mod62_pinvb(uint64_t p)
{
  assert(mod62_valid(p));

  uint64_t pinvb;
  pinvb = (0xf050309070d0b010ULL >> (4 * (p & 15)));   // correct mod 2^4
  pinvb *= (2 - p * pinvb);    // correct mod 2^8
  pinvb *= (2 - p * pinvb);    // correct mod 2^16
  pinvb *= (2 - p * pinvb);    // correct mod 2^32
  pinvb *= (2 - p * pinvb);    // correct mod 2^64
  return pinvb;
}



// in all routines below, pinv and pinvb refer to the values computed by
// mod62_pinv(p) and mod62_pinvb(p)


// c0 in [0, B), c1 in [0, p)
// returns (B*c1 + c0) mod p, in [0, 4p)
static inline uint64_t
mod62_reduce_pinv_lazy4(uint64_t c1, uint64_t c0, uint64_t p, uint64_t pinv)
{
  /*
    Let c = B*c1 + c0, u = floor(4*c/B), q = floor(pinv*u/B) + u.
    Note that q does not overflow because
      p*q <= p*(pinv + B)*u/B < (B^2/4)*(u/B) <= B/4*(4*c/B) <= c <= p*B.
    We have
      0 <= B^2/4p - (pinv + B) < 1,
      0 <= 4*c/B - u < 1,
      0 <= pinv*u/B + u - q < 1.
    Multiplying respectively by p*u/B, B/4, and p, and adding, yields
      0 <= c - q*p < 4p.
  */
  uint64_t u = (c1 << 2) + (c0 >> 62);
  uint64_t q = MUL_HI(pinv, u) + u;
  return c0 - q * p;
}


// same as mod62_reduce_pinv_lazy4(), but returns result in [0, 2p)
static inline uint64_t
mod62_reduce_pinv_lazy2(uint64_t c1, uint64_t c0, uint64_t p, uint64_t pinv)
{
  uint64_t r = mod62_reduce_pinv_lazy4(c1, c0, p, pinv);
  return r - ((r >= 2*p) ? (2*p) : 0);
}


// same as mod62_reduce_pinv_lazy4(), but returns result in [0, p)
static inline uint64_t
mod62_reduce_pinv(uint64_t c1, uint64_t c0, uint64_t p, uint64_t pinv)
{
  uint64_t r = mod62_reduce_pinv_lazy2(c1, c0, p, pinv);
  return r - ((r >= p) ? p : 0);
}


// x, y in [0, B) with x * y < p * B
// returns x * y mod p, in [0, 4p)
static inline uint64_t
mod62_mul_pinv_lazy4(uint64_t x, uint64_t y, uint64_t p, uint64_t pinv)
{
  uint64_t c0, c1;
  MUL_WIDE(c1, c0, x, y);
  return mod62_reduce_pinv_lazy4(c1, c0, p, pinv);
}


// same as mod62_mul_pinv_lazy4(), but returns result in [0, 2p)
static inline uint64_t
mod62_mul_pinv_lazy2(uint64_t x, uint64_t y, uint64_t p, uint64_t pinv)
{
  uint64_t c0, c1;
  MUL_WIDE(c1, c0, x, y);
  return mod62_reduce_pinv_lazy2(c1, c0, p, pinv);
}


// same as mod62_mul_pinv_lazy4(), but returns result in [0, p)
static inline uint64_t
mod62_mul_pinv(uint64_t x, uint64_t y, uint64_t p, uint64_t pinv)
{
  uint64_t c0, c1;
  MUL_WIDE(c1, c0, x, y);
  return mod62_reduce_pinv(c1, c0, p, pinv);
}


// y in [0, p)
// returns floor(y * B / p)
static inline uint64_t mod62_ypinv(uint64_t y, uint64_t p, uint64_t pinv)
{
  /*
    we have
      0 <= B^2/4p - (pinv + B) < 1,
    so
      0 <= B*y/p - pinv*4*y/B - 4*y < 4y/B < 1.
    Let u = floor(4*y*pinv/B) + 4*y. Then u < B, and
      0 <= 4*y*pinv/B + 4*y - u < 1.
    Adding:
      0 <= B*y/p - u < 2.
    Let z = floor(B*y/p). Then
      -1 < z - B*y/p <= 0.
    Adding:
      -1 < z - u < 2,   i.e.  0 <= z - u <= 1.
    So the correct value for z is either u or u + 1.
  */
  assert(y < p);
  uint64_t t = 4 * y;
  uint64_t u = MUL_HI(t, pinv) + t;
  return u + ((-u * p) >= p);    // correction if u is one too small
}



// x in [0, B), y in [0, p), ypinv = floor(y * B / p)
// returns x * y mod p, in [0, 2p)
// (This is Shoup's algorithm from NTL, but skipping the final reduction.)
static inline uint64_t
mod62_mul_ypinv_lazy(uint64_t x, uint64_t y, uint64_t ypinv, uint64_t p)
{
  /*
    we have
      0 <= B*y/p - ypinv < 1,
    so
      0 <= x*y - ypinv*x*p/B < x*p/B < p.
    Let q = floor(x * ypinv / B). Then
      0 <= x*ypinv/B - q < 1
    so
      0 <= x*ypinv*p/B - q*p < p.
    Adding:
      0 <= x*y - q*p < 2p.
  */
  uint64_t q = MUL_HI(x, ypinv);
  return x * y - q * p;
}


// same as mod62_mul_ypinv_lazy(), but with output in [0, p)
static inline uint64_t
mod62_mul_ypinv(uint64_t x, uint64_t y, uint64_t ypinv, uint64_t p)
{
  uint64_t r = mod62_mul_ypinv_lazy(x, y, ypinv, p);
  return r - ((r >= p) ? p : 0);
}


// x, y in [0, B) with x * y < p * B
// returns x * y / B mod p, in [0, 2p)
static inline uint64_t
mod62_mul_pinvb_lazy(uint64_t x, uint64_t y, uint64_t p, uint64_t pinvb)
{
  uint64_t c0, c1;
  MUL_WIDE(c1, c0, x, y);
  return c1 - MUL_HI(c0 * pinvb, p) + p;
}


// same as mod62_mul_pinvb_lazy(), but returns result in [0, p)
static inline uint64_t
mod62_mul_pinvb(uint64_t x, uint64_t y, uint64_t p, uint64_t pinvb)
{
  uint64_t c0, c1;
  MUL_WIDE(c1, c0, x, y);
  int64_t r = c1 - MUL_HI(c0 * pinvb, p);
  return r + ((r < 0) ? p : 0);
}


// x in [0, p), n >= 0
// returns x^n mod p, in [0, p)
// (uses stupid hardware division)
uint64_t mod62_pow(uint64_t x, uint64_t n, uint64_t p);


// x in [0, p), n >= 0
// returns x^n mod p, in [0, p)
uint64_t mod62_pow_pinv(uint64_t x, uint64_t n, uint64_t p, uint64_t pinv);


// returns 2^k mod p in [0, p)
// (negative k is allowed)
uint64_t mod62_2exp(int k, uint64_t p, uint64_t pinv);


// 0 < x, y < 2^64
// returns d = gcd(x, y) > 0, and s, t such that
//    s*x - t*y = d
//    0 <  s <= y/d
//    0 <= t <  x/d
void mod62_xgcd(uint64_t* d, uint64_t* s, uint64_t* t, uint64_t x, uint64_t y);


// n > 0, any x
// returns inverse of x mod n, or 0 if no inverse exists
uint64_t mod62_inv(uint64_t x, uint64_t n);


#endif
