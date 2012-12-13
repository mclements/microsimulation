/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_gsl.h                                                        *
 *                                                                           *
 *   PURPOSE:                                                                *
 *     Function prototypes for using URNG of type GSL:                       *
 *     uniform random number from GSL (GNU Scientific Library),              *
 *     see http://www.gnu.org/software/gsl/.                                 *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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
#ifndef URNG_GSL_H_SEEN
#define URNG_GSL_H_SEEN
/*---------------------------------------------------------------------------*/
#include <gsl/gsl_rng.h>
/*---------------------------------------------------------------------------*/

/* 
   =NODE  URNG-GSL  Interface to GSL uniform random number generators

   =UP URNG [20]

   =DESCRIPTION
      Interface to the uniform random number generators from the
      GNU Scientific Library (GSL). Documentation and source code 
      of this library is available from
      @uref{http://www.gnu.org/software/gsl/}.

      The interface to the GSL must be compiled into UNU.RAN using the
      configure flag @code{--with-urng-gsl}.
      Notice that the GSL has to be installed before running
      @code{./configure}.

   =HOWTOUSE
      When using this interface @file{unuran_urng_gsl.h} must be included
      in the corresponding C file, i.e., one must add the line
      @smallexample
      #include <unuran_urng_gsl.h>
      @end smallexample

      Moreover, one must not forget to link the executable against
      @file{libgsl}.

      The following routines are supported for URNG objects of
      type GSL:

      @itemize @minus
      @item unur_urng_sample()
      @item unur_urng_sample_array()
      @item unur_urng_seed() 
      @item unur_urng_reset() 
      @item unur_urng_free()
      @end itemize

      @smallexample
      @include ref_example_gsl.texi
      @end smallexample

   =END

*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

/*---------------------------------------------------------------------------*/

UNUR_URNG *unur_urng_gsl_new( const gsl_rng_type *urngtype );
/*
   Make object for URNGs from the @file{GSL} (GNU Scientific Library).
   @var{urngtype} is the type of the chosen generator as described in the
   GSL manual (see Section Random Number Generation). This library is
   available from @uref{http://www.gnu.org/software/gsl/}.
*/

UNUR_URNG *unur_urng_gslptr_new( gsl_rng *urng );
/*
   Similar to unur_urng_gsl_new() but it uses a pointer to a
   generator object as returned by @code{gsl_rng_alloc(rng_type)};
   see @file{GSL} manual for details.

   @emph{Notice}: There is a subtle but important difference between
   these two calls. When a generator object is created by a 
   unur_urng_gsl_new() call, then resetting of the generator works.
   When a generator object is created by a unur_urng_gslptr_new()
   call, then resetting only works after a
   @code{unur_urng_seed(urng,myseed)} call. 
*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* URNG_GSL_H_SEEN */
/*---------------------------------------------------------------------------*/
