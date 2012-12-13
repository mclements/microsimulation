/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_randomshift.c                                                *
 *                                                                           *
 *   Meta uniform random number generator for Randomized Quasi-Monte Carlo   *
 *   integration.                                                            *
 *                                                                           * 
 *   (1) Sample and store a random vector S.                                 * 
 *   (2) Run a QMC simulation where S is added to each point of the          *
 *       quasi-random point (mod 1).                                         *
 *   (3) Repeat steps (1) and (2).                                           *
 *   (4) Return sample mean and sample variance.                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   References:                                                             *
 *                                                                           *
 *   [1] Cranley, R. and Patterson T.N.L. (1976),                            *
 *       Randomization of number theoretic methods for multiple integration. *
 *       SIAM J. Num. Anal. 13, 904-914.                                     *
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
#include "urng_randomshift.h"
/*---------------------------------------------------------------------------*/

/* structure for storing generators  */
struct unur_urng_randomshift {
  UNUR_URNG *qrng;   /* generator for point set */
  UNUR_URNG *srng;   /* generator for shift vector */
  double *shift;     /* shift vector */
  double *X;         /* working array for storing points */
  int dim;           /* dimension for QRNG         */
  int n;             /* current coordinate         */
};

/*---------------------------------------------------------------------------*/

static double _unur_urng_randomshift_sample( struct unur_urng_randomshift *rs );
/*---------------------------------------------------------------------------*/
/* Sample next coordinate of random point.                                   */
/* Skip to first coordinate of next point after last coordinate.             */
/*---------------------------------------------------------------------------*/

static int _unur_urng_randomshift_sample_array( struct unur_urng_randomshift *rs, double *X, int dim );
/*---------------------------------------------------------------------------*/
/* Sample random point and store in X. use only the first dim coordinates    */
/* if dim is less than the dimension of the generated points.                */
/*---------------------------------------------------------------------------*/

static void _unur_urng_randomshift_free( struct unur_urng_randomshift *rs );
/*---------------------------------------------------------------------------*/
/* Free generator object.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_urng_randomshift_reset( struct unur_urng_randomshift *rs );
/*---------------------------------------------------------------------------*/
/* Reset generator object.                                                   */
/*---------------------------------------------------------------------------*/

static void _unur_urng_randomshift_nextpoint( struct unur_urng_randomshift *rs );
/*---------------------------------------------------------------------------*/
/* Skip to first coordinate of next point.                                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

double
_unur_urng_randomshift_sample( struct unur_urng_randomshift *rs )
     /*----------------------------------------------------------------------*/
     /* Sample next coordinate of random point.                              */
     /* Skip to first coordinate of next point after last coordinate.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rs  ... pointer to URNG object                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   next coordinate                                                    */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  int i;
  int n = rs->n;
  double *X = rs->X;

  if (n==0) {
    /* first coordinate --> get new point */
    unur_urng_sample_array(rs->qrng,X,rs->dim);
    /* make shift */
    for (i=0; i<rs->dim; i++) {
      X[i] += rs->shift[i];
      if (X[i] >= 1.) X[i] -= 1.;
      if (X[i] < 0.) X[i] += 1.;
    }
  }

  /* increment counter for coordinate */
  rs->n = (n+1) % rs->dim;

  /* return result */
  return X[n];

} /* end of _unur_urng_randomshift_sample() */

/*---------------------------------------------------------------------------*/

int
_unur_urng_randomshift_sample_array( struct unur_urng_randomshift *rs, double *X, int dim )
     /*----------------------------------------------------------------------*/
     /* Sample random point and store in X. use only the first dim coordin.  */
     /* if dim is less than the dimension of the generated points.           */
     /* Return number of entries filled into array X.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rs     ... pointer to URNG object                                  */
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
  int i;

  if (dim > rs->dim) dim = rs->dim;

  if (dim == rs->dim) {
    unur_urng_sample_array(rs->qrng,X,dim); }
  else {  /* dim < rs->dim */
    /* size of array X to small --> use working array */
    unur_urng_sample_array(rs->qrng,rs->X,dim);
    memcpy( X, rs->X, dim*sizeof(double) );
  }

  /* make shift */
  for (i=0; i<dim; i++) {
    X[i] += rs->shift[i];
    if (X[i] >= 1.) X[i] -= 1.;
    if (X[i] < 0.) X[i] += 1.;
  }

  return dim;
} /* end of _unur_urng_randomshift_sample_array() */

/*---------------------------------------------------------------------------*/

void
_unur_urng_randomshift_free( struct unur_urng_randomshift *rs )
     /*----------------------------------------------------------------------*/
     /* Free URNG object.                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rs   ... pointer to URNG object                                    */
     /*----------------------------------------------------------------------*/
{
  if (rs) {
    free( rs->shift );
    free( rs->X );
    free( rs );
  }
} /* end of _unur_urng_randomshift_free() */

/*---------------------------------------------------------------------------*/

void
_unur_urng_randomshift_reset( struct unur_urng_randomshift *rs )
     /*----------------------------------------------------------------------*/
     /* Reset URNG object.                                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rs   ... pointer to URNG object                                    */
     /*----------------------------------------------------------------------*/
{
  unur_urng_reset( rs->qrng );
  unur_urng_reset( rs->srng );
  rs->n = 0u;
} /* end of _unur_urng_randomshift_reset() */

/*---------------------------------------------------------------------------*/

void
_unur_urng_randomshift_nextpoint( struct unur_urng_randomshift *rs )
     /*----------------------------------------------------------------------*/
     /* Skip to first coordinate of next point.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rs   ... pointer to URNG object                                    */
     /*----------------------------------------------------------------------*/
{
  rs->n = 0u;
} /* end of _unur_urng_randomshift_nextpoint() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_randomshift_new( UNUR_URNG *qrng, UNUR_URNG *srng, int dim )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type GSL-QRNG                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   qrng ... generator for point set                                   */
     /*   srng ... generator for random shifts                               */
     /*   dim  ... dimension of point sets                                   */
     /*----------------------------------------------------------------------*/
{
  struct unur_urng_randomshift *rs;
  UNUR_URNG *urng;

  /* check argument */
  _unur_check_NULL( "URNG", qrng, NULL );  COOKIE_CHECK(qrng,CK_URNG,NULL);
  _unur_check_NULL( "URNG", srng, NULL );  COOKIE_CHECK(srng,CK_URNG,NULL);

  /* make structure to store random shift object */
  rs = _unur_xmalloc( sizeof(struct unur_urng_randomshift) );
  rs->shift = _unur_xmalloc( dim * sizeof(double) );
  rs->X = _unur_xmalloc( dim * sizeof(double) );
  rs->qrng = qrng;
  rs->srng = srng;
  rs->dim = dim;
  rs->n = 0u;

  /* make UNURAN_URNG object */
  urng = unur_urng_new ( (double(*)(void*)) _unur_urng_randomshift_sample, rs );
  unur_urng_set_sample_array (urng, (unsigned int(*)(void*,double*,int)) _unur_urng_randomshift_sample_array);
  unur_urng_set_delete (urng, (void(*)(void*)) _unur_urng_randomshift_free);
  unur_urng_set_reset (urng, (void(*)(void*)) _unur_urng_randomshift_reset);
  unur_urng_set_sync (urng, (void(*)(void*)) _unur_urng_randomshift_nextpoint);

  /* initialize random shift */
  unur_urng_sample_array(rs->srng,rs->shift,rs->dim);

  return urng;

} /* end of unur_urng_randomshift_new() */

/*---------------------------------------------------------------------------*/

int
unur_urng_randomshift_nextshift( UNUR_URNG *urng )
     /*----------------------------------------------------------------------*/
     /* Get the next (randomly chosen) vector for shifting the points set.   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rs   ... pointer to URNG object                                    */
     /*----------------------------------------------------------------------*/
{
  struct unur_urng_randomshift *rs = urng->state;
  unur_urng_sample_array(rs->srng,rs->shift,rs->dim);
  unur_urng_reset( rs->qrng );
  rs->n = 0u;
  return UNUR_SUCCESS;
} /* end of unur_urng_randomshift_nextshift() */

/*---------------------------------------------------------------------------*/
#endif   /* #if defined(UNUR_URNG_UNURAN) */
/*---------------------------------------------------------------------------*/
