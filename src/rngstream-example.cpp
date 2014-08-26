
#include "rngstream-boost.hpp"
#include <boost/random/uniform_01.hpp>
#include <iostream>

#include "RngStream.cpp"

int main() {
  boost::rngstream gen[2];
  boost::uniform_01<> U01;
  std::cout << "Expected: 0.127011; observed: " << U01(gen[0]) << std::endl;
  std::cout << "Expected: 0.759582; observed: " << U01(gen[1]) << std::endl;
  gen[0].ResetNextSubstream();
  std::cout << "Expected: 0.079399; observed: " << U01(gen[0]) << std::endl;
  return 0;
}
// g++ -I. -I/home/marcle/R/x86_64-pc-linux-gnu-library/3.1/BH/include rngstream-example.cpp
// R -q -e "require(parallel); base=c(407L,rep(12345L,6)); .Random.seed=base; runif(2); .Random.seed=nextRNGStream(base); runif(2); .Random.seed=nextRNGSubStream(base); runif(2)"
