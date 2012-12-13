/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_fvoid.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *     Function prototypes for using uniform of type FVOID                   *
 *     (i.e. routine without an argment:   double uniform(void)              *
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
#ifndef URNG_FVOID_H_SEEN
#define URNG_FVOID_H_SEEN
/*---------------------------------------------------------------------------*/

/* 
   =NODE  URNG-FVOID  Simple interface for uniform random number generators

   =UP URNG [10]

   =DESCRIPTION
      Simple interface for URNGs of type @code{double uniform(void *state)}.

      UNU.RAN contains some build-in URNGs of this type:
      @table @code
      @item unur_urng_MRG31k3p
      Combined multiple recursive generator by Pierre L'Ecuyer and
      Renee Touzin.
      @item unur_urng_fish
      Linear congruential generator by Fishman and Moore.
      @item unur_urng_mstd
      Linear congruential generator "Minimal Standard" by Park and Miller.
      @end table
      
      Notice, however, that these generators are provided as a
      fallback for the case that no state-of-the-art uniform random
      number generators (e.g. @pxref{URNG-RNGSTREAM,Pierre L'Ecuyer's
      @file{Rngstream} library, Pierre L'Ecuyer's @file{Rngstream}
      library}) are used.
      
   =HOWTOUSE
      Create an URNG object using unur_urng_fvoid_new(). 
      By this call a pointer to the sampling routine and (optional) a
      pointer to a reset routine are copied into the URNG object.
      Other functions, like seeding the URNG, switching to antithetic
      random number, or jumping to next substream, can be added to the
      URNG object by the respective calls, e.g. by
      unur_urng_set_seed().

      The following routines are supported for URNG objects of this
      type: 

      @itemize @minus
      @item unur_urng_sample()
      @item unur_urng_sample_array()
      @item unur_urng_seed()   [optional]
      @item unur_urng_reset()   [optional]
      @item unur_urng_free()
      @end itemize

   =END
*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

UNUR_URNG *unur_urng_fvoid_new( double (*urand)(void *state), void (*reset)(void *state) );
/*
   Make a URNG object for a generator that consists of a single
   function call @var{urand}.

   If there is no @var{reset} function use NULL for the second argument.
*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* URNG_FVOID_H_SEEN */
/*---------------------------------------------------------------------------*/
