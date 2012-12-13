/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_builtin.c                                                    *
 *                                                                           *
 *      routines for default built-in uniform random number generators.      *
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
#include <urng/urng.h>
#include "urng_fvoid.h"
#include "urng_builtin.h"
/*---------------------------------------------------------------------------*/

UNUR_URNG *unur_urng_builtin( void )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of FVOID.                                        */
     /*                                                                      */
     /* parameters: none                                                     */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng;
  urng = unur_urng_fvoid_new(unur_urng_MRG31k3p, unur_urng_MRG31k3p_reset);
  unur_urng_set_seed(urng, unur_urng_MRG31k3p_seed);
  return urng;
} /* unur_urng_builtin() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *unur_urng_builtin_aux( void )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of FVOID.                                        */
     /*                                                                      */
     /* parameters: none                                                     */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng;
  urng = unur_urng_fvoid_new(unur_urng_fish, unur_urng_fish_reset);
  unur_urng_set_seed(urng, unur_urng_fish_seed);
  return urng;
} /* end of unur_urng_builtin_aux() */

/*---------------------------------------------------------------------------*/
