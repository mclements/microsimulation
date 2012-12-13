/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_gsl.c                                                        *
 *                                                                           *
 *   routines to get new URNG object with sampling routine of type GSL.      *
 *   GSL (GNU Scientific Library), see http://www.gnu.org/software/gsl/.     *
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
#if defined(UNURAN_HAS_GSL) && defined(UNUR_URNG_UNURAN)
/*---------------------------------------------------------------------------*/
#include <urng/urng.h>
#include "urng_gsl.h"
/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_gslptr_new( gsl_rng *gsl )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type GSL.                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gsl ... pointer to generator structure                             */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng;

  /* check argument */
  if (gsl == NULL) {
    _unur_error("URNG",UNUR_ERR_NULL,"Cannot create GSL object");
    return NULL;
  }

  urng = unur_urng_new( (double(*)(void*)) gsl_rng_uniform_pos, gsl );
  unur_urng_set_delete(urng, (void(*)(void*)) gsl_rng_free);
  unur_urng_set_seed(urng, (void(*)(void*,unsigned long)) gsl_rng_set);
  return urng;
} /* end of unur_urng_gslptr_new() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_gsl_new( const gsl_rng_type *urngtype )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type GSL.                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urngtype ... type of generator                                     */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = unur_urng_gslptr_new(gsl_rng_alloc(urngtype));
  unur_urng_seed(urng,gsl_rng_default_seed);
  return urng;
} /* end of unur_urng_gsl_new() */

/*---------------------------------------------------------------------------*/
#endif   /* defined(UNURAN_HAS_GSL) && defined(UNUR_URNG_UNURAN) */
/*---------------------------------------------------------------------------*/
