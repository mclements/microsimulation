
#ifndef RNG_WRAPPER_H
#define RNG_WRAPPER_H

extern "C" {

#include "RngStream.h"

  class Rng {
  public:
    Rng(std::string n = "");
    ~Rng();
    void set();
    void nextSubstream();
  private:  
    RngStream stream;
  };

  void r_create_current_stream();

  void r_remove_current_stream();

  void test_rstream2(double * x);

} // extern "C"

#endif
