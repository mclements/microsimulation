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
 *         Linear congruential generator with                                *
 *         m = 2^31-1, a = 742938285, c = 0.                                 *
 *                                                                           *
 *   WARNING:                                                                *
 *         Not state-of-the-art. SHOULD NOT BE USED ANY MORE.                *
 *         In UNU.RAN only as auxilliary second stream.                      *
 *         Should be replaced in future releases.                            *
 *                                                                           *
 *   REFERENCE:                                                              *
 *   Fishman G.S., Moore L.R. (1986): An exhaustive analysis of              *
 *      multiplicative congruential number generators with modulus 2^31-1,   *
 *      SIAM Journal Sci. Stat. Comput., 24-45.                              *
 *                                                                           *
 *   Copyright for generator code by Ernst Stadlober.                        * 
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
#define SEED  (12345L)

/* status variable */
static unsigned long x = SEED;
static unsigned long x_start = SEED; /* seed of last stream */

/*---------------------------------------------------------------------------*/

double
unur_urng_fish (void *dummy ATTRIBUTE__UNUSED)
{
  
# define A   742938285
# define AHI (A>>15)
# define ALO (A&0x7FFF)

  unsigned long xhi, xlo, mid;   /* for intermediate results */

  /* generator */
  xhi = x>>16;
  xlo = x&0xFFFF;
  mid = AHI*xlo + (ALO<<1)*xhi;
  x   = AHI*xhi + (mid>>16) + ALO*xlo;
  if (x&0x80000000) x -= 0x7FFFFFFF;
  x += ((mid&0xFFFF)<<15);
  if (x&0x80000000) x -= 0x7FFFFFFF;

  return (x*4.656612875245797e-10);
} /* end of unur_urng_fish() */

/*---------------------------------------------------------------------------*/

void
unur_urng_fish_seed (void *dummy ATTRIBUTE__UNUSED, unsigned long seed)
{
  if (seed==0) {
    _unur_error("URNG.fish",UNUR_ERR_GENERIC,"seed = 0");
    return;
  }
  
  x_start = seed;
  x = seed;

} /* end of unur_urng_fish_seed() */

/*---------------------------------------------------------------------------*/

void
unur_urng_fish_reset (void *dummy ATTRIBUTE__UNUSED)
{
  x = x_start;
} /* end of unur_urng_fish_reset() */

/*---------------------------------------------------------------------------*/
