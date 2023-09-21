/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/


unsigned zz_gcd(unsigned x, unsigned y)
{
  // couple of special cases, since we often call with y == num_primes,
  // which is a small positive integer
  switch (y)
    {
    case 0: return x;
    case 1: return 1;
    case 2: return (x % 2) ? 1 : 2;
    case 3: return (x % 3) ? 1 : 3;
    case 4:
      switch (x % 4)
	{
	case 0: return 4;
	case 1: return 1;
	case 2: return 2;
	case 3: return 1;
	}
    default:
      return zz_gcd(y, x % y);
    }
}
