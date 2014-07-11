#include "microsimulation.hpp"

namespace ssim {

  extern "C" {
    
  /** 
      @brief A utility function to create the current_stream.
      Used when initialising the microsimulation package in R.
  */
    void r_create_current_stream()
    {
      //*default_stream() = *(new Rng());
      //*current_stream() = *default_stream();
    }
    
  /** 
      @brief A utility function to remove the current_stream.
      Used when finalising the microsimulation package in R.
  */
    void r_remove_current_stream()
    {
      //delete default_stream();
    }
    
  /** 
      @brief A utility function to set the user random seed for the simulation.
  */
    void r_set_user_random_seed(double * inseed) {
      unsigned long seed[6];
      for(int i=0; i<6; i++) {
  	seed[i] = (unsigned long)inseed[i];
      }
      Rng::SetPackageSeed(seed);
      default_stream()->SetSeed(seed);
    }
    
  /** 
      @brief A utility function to set the user random seed for the simulation.
  */
    void r_get_user_random_seed(double * outseed) {
      unsigned long seed[6];
      default_stream()->GetState(seed);
      for(int i=0; i<6; i++) {
  	outseed[i] = (double)seed[i];
      }
    }
    
  /** 
      @brief A utility function to move to the next user random stream for the simulation.
  */
    void r_next_rng_substream() {
      default_stream()->ResetNextSubstream(); // From R, this assumes that default_stream == current_stream.
    }
    
    double *user_unif_rand ()
    {
      if (!current_stream()) {
  	REprintf("user_unif_rand(): No stream created yet!");
  	return NULL;
      }
      *rn() = current_stream()->RandU01();
      return rn();
    }
    
  /** 
      @brief Simple test of the random streams (with a stupid name)
  */
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
