/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_uniform.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *      function prototypes for built-in uniform random number generators    *
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
#ifndef URNG_BUILTIN_H_SEEN
#define URNG_BUILTIN_H_SEEN
/*---------------------------------------------------------------------------*/

/* Combined multiple recursive generator by Pierre L'Ecuyer and Renee Touzin */
/* Copyright (c) 2002 Renee Touzin.                                          */

double unur_urng_MRG31k3p (void *dummy);
void unur_urng_MRG31k3p_seed (void *dummy, unsigned long seed);
void unur_urng_MRG31k3p_reset (void *dummy);

/* Linear congruential generator by Fishman and Moore                        */
/* m = 2^31-1, a = 742938285, c = 0.                                         */

double unur_urng_fish (void *dummy);
void unur_urng_fish_seed (void *dummy, unsigned long seed);
void unur_urng_fish_reset (void *dummy);

/* Linear congruential generator "Minimal Standard"                          */
/* m = 2^31-1, a = 16807, c = 0.                                             */

double unur_urng_mstd (void *dummy);
void unur_urng_mstd_seed (void *dummy, unsigned long seed);
void unur_urng_mstd_reset (void *dummy);

/*---------------------------------------------------------------------------*/

UNUR_URNG *unur_urng_builtin( void );
/*
   Make object for the default builtin URNG.
*/

UNUR_URNG *unur_urng_builtin_aux( void );
/*
   Make object for the default builtin URNG.
*/

/*---------------------------------------------------------------------------*/
#endif  /* URNG_BUILTIN_H_SEEN */
/*---------------------------------------------------------------------------*/
