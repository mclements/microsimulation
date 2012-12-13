/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_fvoid.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *     Function prototypes for using uniform of type RNGSTREAMSPRNG:         *
 *     Pierre L'Ecuyer's RNGSTREAMS package.                                 *
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
#ifndef URNG_RNGSTREAMS_H_SEEN
#define URNG_RNGSTREAMS_H_SEEN
/*---------------------------------------------------------------------------*/
#include <RngStream.h>
/*---------------------------------------------------------------------------*/

/* 
   =NODE  URNG-RNGSTREAM  Interface to L'Ecuyer's RNGSTREAM random number generators

   =UP URNG [50]

   =DESCRIPTION
      URNGs from Pierre L'Ecuyer's @file{RngStream} library for multiple 
      independent streams of pseudo-random numbers. 
      This library provides multiple independent streams of
      pseudo-random numbers which itselves can be splitted into many
      substreams. It is available from
      @uref{http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/}.
      A GNU-style package is available from
      @uref{http://statmath.wu.ac.at/software/RngStreams/}.

      The interface to the RngStream library must be compiled into UNU.RAN using the
      configure flag @code{--with-urng-rngstream}.
      Notice that the RngStream library has to be installed before running
      @code{./configure}.

   =HOWTOUSE
      When using this interface @file{unuran_urng_rngstream.h} must be included
      in the corresponding C file, i.e., one must add the line
      @smallexample
      #include <unuran_urng_rngstream.h>
      @end smallexample

      Moreover, one must not forget to link the executable against the
      @file{RngStream} library (i.e., when using the GNU-style package
      in UNIX like environments one has to add @code{-lrngstreams}
      when linking an executable).

      Notice that the @file{rngstream} library uses a package seed,
      that means one should seed the uniform random number generator
      only once in an application using the routine
      @uref{http://statmath.wu.ac.at/software/RngStreams/doc/rngstreams.html#index-RngStream_005fSetPackageSeed-2,,@code{RngStream_SetPackageSeed}}:

      @smallexample
        unsigned long seed[] = @{111u, 222u, 333u, 444u, 555u, 666u@};
        RngStream_SetPackageSeed(seed);
      @end smallexample

      The following routines are supported for URNG objects of this
      type:

      @itemize @minus
      @item unur_urng_sample()
      @item unur_urng_sample_array()
      @item unur_urng_reset() 
      @item unur_urng_nextsub() 
      @item unur_urng_resetsub() 
      @item unur_urng_anti() 
      @item unur_urng_free()
      @end itemize

   =END

*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

/*---------------------------------------------------------------------------*/

UNUR_URNG *unur_urng_rngstream_new( const char *urngstr );
/*
   Make object for URNGs from Pierre L'Ecuyer's @file{RngStream}
   library. @var{urngstr} is an arbitrary string to label a stream. It
   need not be unique.
*/

UNUR_URNG *unur_urng_rngstreamptr_new( RngStream rngstream );
/*
   Similar to unur_urng_rngstream_new() but it uses a pointer to a 
   generator object as returned by @code{RngStream_CreateStream()}.
*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* URNG_RNGSTREAMS_H_SEEN */
/*---------------------------------------------------------------------------*/
