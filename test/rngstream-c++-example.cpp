
#include "rngstream-c++11.hpp"
#include <random>
#include <iostream>

int main() {
  rngstream gen, gen2;
  std::uniform_real_distribution<> dist(0.0,1.0);
  std::cout << gen2 << std::endl;
  std::cout << "Expected: 0.127011; observed: " << dist(gen) << std::endl;
  std::cout << "Expected: 0.759582; observed: " << dist(gen2) << std::endl;
  gen.ResetNextSubstream();
  std::cout << "Expected: 0.079399; observed: " << dist(gen) << std::endl;
  std::cout << gen << std::endl;
  return 0;
}
// R -q -e "require(parallel); base=c(407L,rep(12345L,6)); .Random.seed=base; runif(1); .Random.seed=nextRNGStream(base); runif(1); .Random.seed=nextRNGSubStream(base); runif(1)"
