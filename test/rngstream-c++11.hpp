#ifndef CPP_RANDOM_RNGSTREAM_HPP
#define CPP_RANDOM_RNGSTREAM_HPP

#include <iostream>
#include <stdexcept>
#include <random>
#include <string>
#include <cstdlib>

#include "RngStream-revised.h"
#include "RngStream-revised.cpp"

class rngstream : public RngStream
{
public:
  typedef unsigned long result_type; 

  /**
   * Returns the smallest value that the generator can produce
   */
  static constexpr result_type min () { return 0ul; }
  /**
   * Returns the largest value that the generator can produce
   */
  static constexpr result_type max ()
  { return 4294967088ul; } // m1+1
  
  /** Seeds the generator with the default seed. */
  rngstream() : RngStream("") { }
  
  // compiler-generated copy ctor and assignment operator are fine
  
  /**  Returns the next value of the generator. */
  result_type operator()() { return U(); }
        
  /**  Writes a @c rngstream to a @c std::ostream. */
  template<class CharT,class Traits>
  friend std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT,Traits>& os, const rngstream& r)
  { 
    result_type seed[6];
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
    result_type seed[6];
    for (int i = 0; i<6; i++)
      is >> seed[i]; 
    r.SetSeed(seed);
    return is; 
  }
  
  /**
   * Returns true if the two generators will produce identical
   * sequences of values.
   */
  friend bool operator==(const rngstream& x, const rngstream& y)
  { 
    result_type seedx[6], seedy[6];
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

#endif // CPP_RANDOM_RNGSTREAM_HPP
