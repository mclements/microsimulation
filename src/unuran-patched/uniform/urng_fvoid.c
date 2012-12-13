/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_fvoid.c                                                      *
 *                                                                           *
 *   routines to get new URNG object with sampling routine of type FVOID:    *
 *   double(*urng)(void), NULL) and global state variable.                   *
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
/*---------------------------------------------------------------------------*/
#if defined(UNUR_URNG_UNURAN)
/*---------------------------------------------------------------------------*/
#include <urng/urng.h>
#include "urng_fvoid.h"
/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_fvoid_new( double (*urand)(void *state), void (*reset)(void *state) )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type FVOID                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urand   ... pointer to uniform random number generator             */
     /*   reset   ... pointer to reset function for URNG                     */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = unur_urng_new( urand, NULL );
  unur_urng_set_reset( urng, reset );
  return urng;
} /* end of unur_urng_fvoid_new() */

/*---------------------------------------------------------------------------*/
#endif   /* #if defined(UNUR_URNG_UNURAN) */
/*---------------------------------------------------------------------------*/

