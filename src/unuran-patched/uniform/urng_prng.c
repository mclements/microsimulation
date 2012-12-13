/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_prng.c                                                       *
 *                                                                           *
 *   routines to get new URNG object with sampling routine of type PRNG.     *
 *   (Lendl's prng package, see http://statmath.wu.ac.at/prng/ or            *
 *   http://random.mat.sbg.ac.at/.                                           *
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
#if defined(UNURAN_HAS_PRNG) && defined(UNUR_URNG_UNURAN)
/*---------------------------------------------------------------------------*/
#include <urng/urng.h>
#include "urng_prng.h"
/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_prngptr_new( struct prng *prng )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type PRNG.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   prng ... pointer to generator structure                            */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng;

  if (prng==NULL) {
    _unur_error("URNG",UNUR_ERR_NULL,"Cannot create PRNG object");
    return NULL;
  }

  urng = unur_urng_new( (double(*)(void*)) prng->get_next, prng );
  unur_urng_set_reset(urng, (void(*)(void*)) prng->reset);
  unur_urng_set_delete(urng, (void(*)(void*)) prng->destroy);

  if (prng_can_seed(prng)) {
    unur_urng_set_seed(urng, (void(*)(void*,unsigned long))prng->seed);
  }

  return urng;
} /* end of unur_urng_prngptr_new() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_prng_new( const char *prngstr )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type PRNG.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   prngstr ... string that describes generator                        */
     /*----------------------------------------------------------------------*/
{
  struct prng *prng = prng_new(prngstr);
  if (prng==NULL) {
    _unur_error("URNG",UNUR_ERR_NULL,"Cannot create PRNG object for given string");
    return NULL;
  }
  return unur_urng_prngptr_new (prng);
} /* end of unur_urng_prng_new() */

/*---------------------------------------------------------------------------*/
#endif   /* defined(UNURAN_HAS_PRNG) && defined(UNUR_URNG_UNURAN) */
/*---------------------------------------------------------------------------*/
