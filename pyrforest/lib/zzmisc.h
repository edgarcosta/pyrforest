/*
  Copyright (C) 2013, David Harvey
  See the file COPYING for license details.
*/

#ifndef ZZ_MISC_H
#define ZZ_MISC_H


#include <gmp.h>


#define MIN(xxx, yyy) (((xxx) < (yyy)) ? (xxx) : (yyy))
#define MAX(xxx, yyy) (((xxx) > (yyy)) ? (xxx) : (yyy))


// used for splitting threads into teams, for load balancing
unsigned zz_gcd(unsigned x, unsigned y);


#endif
