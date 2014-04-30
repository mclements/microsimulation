
#include "rngstream-boost.hpp"
#include <boost/random/uniform_01.hpp>
#include <iostream>

int main() {
  boost::rngstream gen, gen2;
  boost::uniform_01<> dist;
  std::cout << gen2 << std::endl;
  std::cout << "Expected: 0.127011; observed: " << dist(gen) << std::endl;
  std::cout << "Expected: 0.759582; observed: " << dist(gen2) << std::endl;
  gen.ResetNextSubstream();
  std::cout << "Expected: 0.079399; observed: " << dist(gen) << std::endl;
  std::cout << gen << std::endl;
  return 0;
}
// R -q -e "require(parallel); base=c(407L,rep(12345L,6)); .Random.seed=base; runif(2); .Random.seed=nextRNGStream(base); runif(2); .Random.seed=nextRNGSubStream(base); runif(2)"
