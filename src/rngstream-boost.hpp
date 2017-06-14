 /* boost random/rngstream-boost.hpp header file
 *
 * Copyright Mark Clements 2014
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org for most recent version including documentation.
 *
 *
 * Revision history
 */

/***********************************************************************\
 *
 * File:           RngStream.cpp for multiple streams of Random Numbers
 * Language:       C++ (ISO 1998)
 * Copyright:      Pierre L'Ecuyer, University of Montreal
 * Notice:         This code can be used freely for personal, academic,
 *                 or non-commercial purposes. For commercial purposes, 
 *                 please contact P. L'Ecuyer at: lecuyer@iro.umontreal.ca
 * Date:           14 August 2001
 *
\***********************************************************************/

#ifndef BOOST_RANDOM_RNGSTREAM_BOOST_HPP
#define BOOST_RANDOM_RNGSTREAM_BOOST_HPP

#include <iostream>
#include <stdexcept>
#include <boost/assert.hpp>
#include <boost/config.hpp>
#include <boost/cstdint.hpp>
#include <boost/limits.hpp>
#include <boost/static_assert.hpp>
#include <boost/integer/static_log2.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/random/detail/config.hpp>
#include <boost/random/detail/const_mod.hpp>
#include <boost/random/detail/seed.hpp>
//#include <boost/random/detail/seed_impl.hpp>
#include <boost/detail/workaround.hpp>
#include <string>
#include <cstdlib>

#include <boost/random/detail/disable_warnings.hpp>

#include <RngStream.h>
#include <RngStream.cpp>

namespace boost {
  namespace random {
    
    using ssim::RngStream;

    class rngstream : public RngStream
    {
    public:
      typedef double result_type;
      
      //BOOST_STATIC_CONSTANT(bool, has_fixed_range = false);
      /**
       * Returns the smallest value that the generator can produce
       */
      static double min BOOST_PREVENT_MACRO_SUBSTITUTION () { return 0.0; }
      /**
       * Returns the largest value that the generator can produce
       */
      static double max BOOST_PREVENT_MACRO_SUBSTITUTION ()
      { return 1.0; } // or 2^24?
      
      /** Seeds the generator with the default seed. */
      rngstream() : RngStream("") { }
      
      // compiler-generated copy ctor and assignment operator are fine
      
      /** Seeds the generator with the default seed. */
      void seed() { 
	unsigned long seed[6];
	for (int i = 0; i<6; i++) 
	  seed[i] = 12345ul;
	SetSeed(seed);
      }
      
      /**  Returns the next value of the generator. */
      double operator()() { return RandU01(); }
      
      /** Fills a range with random values */
      template<class Iter>
      void generate(Iter first, Iter last)
      {
          for(; first != last; ++first) {
              *first = (*this)();
          }
      }
      
#ifndef BOOST_RANDOM_NO_STREAM_OPERATORS
      /**  Writes a @c rngstream to a @c std::ostream. */
      template<class CharT,class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<<(std::basic_ostream<CharT,Traits>& os, const rngstream& r)
      { 
	unsigned long seed[6];
	r.GetState (seed);
	for (int i = 0; i<5; i++)
	  os << seed[i] << ' ';
	os << seed[5];
	return os; 
      }
      
      /** Reads a @c rngstream from a @c std::istream. */
      template<class CharT,class Traits>
      friend std::basic_istream<CharT,Traits>&
      operator>>(std::basic_istream<CharT,Traits>& is, rngstream& r)
      { 
	unsigned long seed[6];
	for (int i = 0; i<6; i++)
	  is >> seed[i]; 
	r.SetSeed(seed);
	return is; 
      }
#endif
      
      /**
       * Returns true if the two generators will produce identical
       * sequences of values.
       */
      friend bool operator==(const rngstream& x, const rngstream& y)
      { 
	unsigned long seedx[6], seedy[6];
	x.GetState (seedx);
	y.GetState (seedy);
	for (int i = 0; i<6; i++)
	  if (seedx[i] != seedy[i])
	    return false;
	return true;
      }
      /**
       * Returns true if the two generators will produce different
       * sequences of values.
       */
      friend bool operator!=(const rngstream& x, const rngstream& y)
      { return !(x == y); }
    };
    
  } // namespace random
  
  using random::rngstream;
  
} // namespace boost

#include <boost/random/detail/enable_warnings.hpp>

#endif // BOOST_RANDOM_RNGSTREAM_BOOST_HPP
