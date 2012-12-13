/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: fish.c                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         uniform random number generator provided by UNU.RAN               *
 *         random number generators inside UNU.RAN.                          *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *         Linear congruential generator                                     *
 *         x_(n+1) = 16807 * x_n mod (2^32 - 1)  ["Minimal Standard"]        *
 *                                                                           *
 *   WARNING:                                                                *
 *         Not state-of-the-art. SHOULD NOT BE USED ANY MORE.                *
 *         In UNU.RAN only as auxilliary second stream.                      *
 *         Should be replaced in future releases.                            *
 *                                                                           *
 *   REFERENCE:                                                              *
 *   Park, S. K. and Miller, K. W. (1988):                                   *
 *      Random number generators: good ones are hard to find,                *
 *      Comm. ACM 31, 1192--1201.                                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#include <unur_source.h>
#include "urng_builtin.h"
/*---------------------------------------------------------------------------*/

/* seed (must not be 0!) */
#define SEED  (1804289L)

/* status variable */
static unsigned long x = SEED;
static unsigned long x_start = SEED; /* seed of last stream */

/*---------------------------------------------------------------------------*/

double
unur_urng_mstd (void *dummy ATTRIBUTE__UNUSED)
{

# define a 16807       /* multiplicator */
# define m 2147483647  /* modulus */
# define q 127773      /* m / a */
# define r 2836        /* m % a */

  int hi, lo, test;   /* intermediate results */

  hi = x / q;
  lo = x % q;
  test = a * lo - r * hi;
  x = (test > 0 ) ? test : test + m;
  return (x * 4.656612875245796924105750827e-10);

} /* end of unur_urng_mstd() */

/*---------------------------------------------------------------------------*/

void
unur_urng_mstd_seed (void *dummy ATTRIBUTE__UNUSED, unsigned long seed)
{
  if (seed==0) {
    _unur_error("URNG.mstd",UNUR_ERR_GENERIC,"seed = 0");
    return;
  }
  
  x = x_start = seed;

} /* end of unur_urng_mstd_seed() */

/*---------------------------------------------------------------------------*/

void
unur_urng_mstd_reset (void *dummy ATTRIBUTE__UNUSED)
{
  x = x_start;
} /* end of unur_urng_mstd_reset() */

/*---------------------------------------------------------------------------*/
