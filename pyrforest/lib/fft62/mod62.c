/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#include "mod62.h"


uint64_t mod62_pow(uint64_t x, uint64_t n, uint64_t p)
{
  // binary powering

  uint64_t acc = 1;       // accumulator
  uint64_t u = x;         // x^(2^i)

  while (n != 0)
    {
      if (n & 1)
	acc = mod62_mul(acc, u, p);
      u = mod62_mul(u, u, p);
      n >>= 1;
    }

  return acc;
}


uint64_t mod62_pow_pinv(uint64_t x, uint64_t n, uint64_t p, uint64_t pinv)
{
  if (n == 0)
    return 1;

  if (n == 1)
    return x;

  if (n == 2)
    return mod62_mul_pinv(x, x, p, pinv);

  uint64_t u = x;         // x^(2^i), in [0, 2p)

  // reduce to case n odd

  for (; !(n & 1); n >>= 1)
    u = mod62_mul_pinv_lazy2(u, u, p, pinv);

  // handle rest of bits of n

  uint64_t acc = u;       // accumulator, in [0, 2p)

  for (n >>= 1; n != 0; n >>= 1)
    {
      u = mod62_mul_pinv_lazy2(u, u, p, pinv);
      if (n & 1)
	acc = mod62_mul_pinv_lazy2(acc, u, p, pinv);
    }

  return mod63_reduce2(acc, p);
}


uint64_t mod62_2exp(int k, uint64_t p, uint64_t pinv)
{
  // todo: surely we can do better than this!!!
  if (k >= 0)
    return mod62_pow_pinv(2, k, p, pinv);
  else
    return mod62_pow_pinv((p + 1) / 2, -k, p, pinv);
}


void mod62_xgcd(uint64_t* d, uint64_t* s, uint64_t* t, uint64_t x, uint64_t y)
{
  uint64_t a = x;
  uint64_t b = y;
  int64_t u00 = 1;
  int64_t u01 = 0;
  int64_t u10 = 0;
  int64_t u11 = 1;
  int sign = 1;

  while (b != 0)
    {
      uint64_t q = a / b;
      uint64_t r = a - q * b;
      a = b;
      b = r;

      int64_t t0 = u00 - q * u10;
      u00 = u10;
      u10 = t0;

      int64_t t1 = u01 - q * u11;
      u01 = u11;
      u11 = t1;

      sign = -sign;
    }

  int64_t _s = u00;
  int64_t _t = -u01;

  if (_t < 0)
    {
      if (sign == 1)
	_s -= u10, _t += u11;
      else
	_s += u10, _t -= u11;
    }

  *d = a;
  *s = _s;
  *t = _t;
}



uint64_t mod62_inv(uint64_t x, uint64_t n)
{
  assert(n > 0);
  if (n == 1 || x == 0)
    return 0;
  uint64_t d, s, t;
  mod62_xgcd(&d, &s, &t, x, n);    // d = s*x - t*y
  return (d > 1) ? 0 : s;
}
