/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_gslqmc.c                                                     *
 *                                                                           *
 *   routines to get new URNG object (of type GSL-QRNG) to use               *
 *   quasi-random sequences from the GSL (GNU Scientific Library),           *
 *   see http://www.gnu.org/software/gsl/.                                   *
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
#include "urng_gslqrng.h"
/*---------------------------------------------------------------------------*/

/* structure for storing QRNG */
struct unur_urng_gslqrng {
  gsl_qrng *qrng;    /* pointer to GSL QRNG obejct */
  double *X;         /* working array for storing points */
  unsigned dim;      /* dimension for QRNG         */
  unsigned n;        /* current coordinate         */ 
};

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* We need wrapper functions for the GSL API.                                */

static UNUR_URNG *_unur_urng_gslqrngptr_new( gsl_qrng *qrngptr, unsigned int dim );
/*---------------------------------------------------------------------------*/
/* Get new URNG object of type GSL-QRNG.                                     */
/*---------------------------------------------------------------------------*/

static double _unur_urng_gslqrng_sample( struct unur_urng_gslqrng *qrng );
/*---------------------------------------------------------------------------*/
/* Sample next coordinate of random point.                                   */
/* Skip to first coordinate of next point after last coordinate.             */
/*---------------------------------------------------------------------------*/

static int _unur_urng_gslqrng_sample_array( struct unur_urng_gslqrng *qrng, double *X, unsigned dim );
/*---------------------------------------------------------------------------*/
/* Sample random point and store in X. use only the first dim coordinates    */
/* if dim is less than the dimension of the generated points.                */
/*---------------------------------------------------------------------------*/

static void _unur_urng_gslqrng_free( struct unur_urng_gslqrng *qrng );
/*---------------------------------------------------------------------------*/
/* Free generator object.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_urng_gslqrng_reset( struct unur_urng_gslqrng *qrng );
/*---------------------------------------------------------------------------*/
/* Reset generator object.                                                   */
/*---------------------------------------------------------------------------*/

static void _unur_urng_gslqrng_nextpoint( struct unur_urng_gslqrng *qrng );
/*---------------------------------------------------------------------------*/
/* Skip to first coordinate of next point.                                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

double
_unur_urng_gslqrng_sample( struct unur_urng_gslqrng *qrng )
     /*----------------------------------------------------------------------*/
     /* Sample next coordinate of random point.                              */
     /* Skip to first coordinate of next point after last coordinate.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   qrng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   next coordinate                                                    */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  int n = qrng->n;

  if (n==0)
    /* first coordinate --> get new point */
    gsl_qrng_get( qrng->qrng, qrng->X );

  /* increment counter for coordinate */
  qrng->n = (n+1) % qrng->dim;

  /* return result */
  return qrng->X[n];
} /* end of _unur_urng_gslqrng_sample() */

/*---------------------------------------------------------------------------*/

int
_unur_urng_gslqrng_sample_array( struct unur_urng_gslqrng *qrng, double *X, unsigned dim )
     /*----------------------------------------------------------------------*/
     /* Sample random point and store in X. use only the first dim coordin.  */
     /* if dim is less than the dimension of the generated points.           */
     /* Return number of entries filled into array X.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   qrng   ... pointer to URNG object                                  */
     /*   X      ... pointer to array of length dim                          */
     /*   dim    ... (maximal) number of entries in X to be set              */
     /*                                                                      */
     /* return:                                                              */
     /*   number of uniform random numbers filled into array                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  if (dim < qrng->dim) {
    /* size of array X to small --> use working array */
    gsl_qrng_get( qrng->qrng, qrng->X );
    memcpy( X, qrng->X, dim*sizeof(double) );
    return dim;
  }
  else {
    gsl_qrng_get( qrng->qrng, X );
    return qrng->dim;
  }
} /* end of _unur_urng_gslqrng_sample_array() */

/*---------------------------------------------------------------------------*/

void
_unur_urng_gslqrng_free( struct unur_urng_gslqrng *qrng )
     /*----------------------------------------------------------------------*/
     /* Free URNG object.                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   qrng   ... pointer to URNG object                                  */
     /*----------------------------------------------------------------------*/
{
  if (qrng) {
    gsl_qrng_free( qrng->qrng );
    free( qrng->X );
    free (qrng);
  }
} /* end of _unur_urng_gslqrng_free() */

/*---------------------------------------------------------------------------*/

void
_unur_urng_gslqrng_reset( struct unur_urng_gslqrng *qrng )
     /*----------------------------------------------------------------------*/
     /* Reset URNG object.                                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   qrng   ... pointer to URNG object                                  */
     /*----------------------------------------------------------------------*/
{
  gsl_qrng_init( qrng->qrng );
  qrng->n = 0u;
} /* end of _unur_urng_gslqrng_reset() */

/*---------------------------------------------------------------------------*/

void
_unur_urng_gslqrng_nextpoint( struct unur_urng_gslqrng *qrng )
     /*----------------------------------------------------------------------*/
     /* Skip to first coordinate of next point.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   qrng   ... pointer to URNG object                                  */
     /*----------------------------------------------------------------------*/
{
  qrng->n = 0u;
} /* end of _unur_urng_gslqrng_nextpoint() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
_unur_urng_gslqrngptr_new( gsl_qrng *qrngptr, unsigned int dim )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type GSL-QRNG                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   qrngptr ... pointer to generator structure                         */
     /*   dim     ... dimension of point sets                                */
     /*----------------------------------------------------------------------*/
{
  struct unur_urng_gslqrng *qrng;
  UNUR_URNG *urng;

  /* check argument */
  if (qrngptr == NULL) {
    _unur_error("URNG",UNUR_ERR_NULL,"Cannot create GSL-QRNG object");
    return NULL;
  }

  /* make structure to store QRNG object */
  qrng = _unur_xmalloc( sizeof(struct unur_urng_gslqrng) );
  qrng->X = _unur_xmalloc( dim * sizeof(double) );
  qrng->qrng = qrngptr;
  qrng->dim = dim;
  qrng->n = 0u;

  /* make UNURAN_URNG object */
  urng = unur_urng_new ( (double(*)(void*)) _unur_urng_gslqrng_sample, qrng );
  unur_urng_set_sample_array (urng, (unsigned int(*)(void*,double*,int)) _unur_urng_gslqrng_sample_array);
  unur_urng_set_delete (urng, (void(*)(void*)) _unur_urng_gslqrng_free);
  unur_urng_set_reset (urng, (void(*)(void*)) _unur_urng_gslqrng_reset);
  unur_urng_set_sync (urng, (void(*)(void*)) _unur_urng_gslqrng_nextpoint);

  return urng;
} /* end of _unur_urng_gslqrngptr_new() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_gslqrng_new( const gsl_qrng_type *qrngtype, unsigned int dim )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type GSL-QRNG                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   qrngtype ... type of generator                                      */
     /*   dim     ... dimension of point sets                                */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (dim == 0u) {
    _unur_error("URNG",UNUR_ERR_GENERIC,"Cannot create GSL-QRNG object for dimension 0");
    return NULL;
  }

  return _unur_urng_gslqrngptr_new( gsl_qrng_alloc(qrngtype, dim), dim );
} /* end of unur_urng_gslqrng_new() */

/*---------------------------------------------------------------------------*/
#endif   /* defined(UNURAN_HAS_GSL) && defined(UNUR_URNG_UNURAN) */
/*---------------------------------------------------------------------------*/
