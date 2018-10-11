/**
 * @file    RngStream.h for multiple streams of Random Numbers
 * @author  Pierre L'Ecuyer, University of Montreal
 * Original date: 14 August 2001
 * Modified by Mark Clements <mark.clements@ki.se> 2014-03-22 for the microsimulation package.
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
*/

 
#ifndef RNGSTREAM_H
#define RNGSTREAM_H
 
#include <string>
#include <stdint.h>

namespace ssim {

class RngStream
{
public:

RngStream (const char *name = "");


static bool SetPackageSeed (const double seed[6]);


void ResetStartStream ();


void ResetStartSubstream ();


void ResetNextSubstream ();


void SetAntithetic (bool a);


void IncreasedPrecis (bool incp);


bool SetSeed (const double seed[6]);


void GenAdvanceState (int32_t e, int32_t c,
		      const double A1[3][3], const double A2[3][3],
		      const double InvA1[3][3], const double InvA2[3][3]);


void AdvanceState (int32_t e, int32_t c);


void AdvanceSubstream (int32_t e, int32_t c);


void AdvanceStream (int32_t e, int32_t c);


void CalcMatrix (int32_t e, int32_t c, double C1[3][3], double C2[3][3]);


void GetState (double seed[6]) const;


/* void WriteState () const; */


/* void WriteStateFull () const; */


double RandU01 ();


int RandInt (int i, int j);



private:

double Cg[6], Bg[6], Ig[6];


bool anti, incPrec;


std::string name;


static double nextSeed[6];


double U01 ();


double U01d ();


};

} // ssim namespace
 
#endif
 

