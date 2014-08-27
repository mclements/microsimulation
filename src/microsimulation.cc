
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

#include "microsimulation.h"

namespace ssim {

  double rweibullHR(double shape, double scale, double hr){
    return R::rweibull(shape, scale*pow(hr,1.0/shape));
  }
  
  static Rng * default_stream, * current_stream;
  static double rn = 0.0;
  
  Rng::~Rng() {
    if (current_stream->id == this->id) // Is this the current stream?
      current_stream = default_stream;  // If so, change the current stream to being the default
  }
  
  void Rng::set() {
    current_stream = this; // make this the current stream
  }

  Rng * Rng::get_default_stream() { return default_stream; }
  Rng * Rng::get_current_stream() { return current_stream; }
  
  extern "C" {
    
    void r_create_current_stream()
    {
      default_stream = new Rng();
      current_stream = default_stream;
    }
    
    void r_remove_current_stream()
    {
      delete default_stream;
    }
    
    void r_set_user_random_seed(double * seed) {
      Rng::SetPackageSeed(seed); // sets the package seed and default stream's seed with the same value 
      default_stream->SetSeed(seed);
    }
    
    void r_get_user_random_seed(double * seed) {
      default_stream->GetState(seed);
    }
    
    void r_next_rng_substream() {
      default_stream->ResetNextSubstream();
    }
    
    double *user_unif_rand ()
    {
      if (!current_stream) {
	REprintf("user_unif_rand(): No stream created yet!");
	return NULL;
      }
      rn = current_stream->RandU01();
      return &rn;
    }
    
    void test_rstream2(double * x) {
      Rng * s1 = new Rng();
      Rng * s2 = new Rng();
      x[0]=WithRNG(s1,R::rexp(1.0));
      x[1]=WithRNG(s2,R::rexp(1.0));
      s1->ResetNextSubstream();
      x[2]=R::rexp(1.0);
      delete s1;
      delete s2;
    }
    
  } // extern "C"
  
} // namespace ssim
  
namespace R {
  double rnormPos(double mean, double sd) {
    double x;
    while ((x=R::rnorm(mean,sd))<0.0) { }
    return x;
  }
}
