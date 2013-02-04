
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Random.h>
//#include "RngStream.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Rmath.h>
#include <string>

extern "C" {

#include "RngStream.h"
#include "Rng_wrapper.h"

#define WithRNG(rng,expr) (rng->set(), expr)

  static RngStream current_stream;
  static double rn = 0.0;

  void r_create_current_stream()
  {
    current_stream = (RngStream) malloc (sizeof (struct RngStream_InfoState));
    
    if (current_stream == NULL) {
      error("r_create_current_stream: No more memory\n\n");
    }
    current_stream->name=(char *)"";
    return;
  }
  
  void r_remove_current_stream()
  {
    free(current_stream);
    /*    RngStream_DeleteStream(&current_stream);*/
    return;
  }

  Rng::Rng(std::string n) {
    stream = RngStream_CreateStream(n.c_str());
  }

  Rng::~Rng() {
    RngStream_DeleteStream(&stream);
  }

  void Rng::set() {
    if (stream == NULL)
      error("B: Invalid NULL pointer");
    current_stream = stream;
    return;
  }

  void Rng::nextSubstream() {
    set(); // is this useful?
    RngStream_ResetNextSubstream(stream);
  }

  double *user_unif_rand ()
  {
    if (!current_stream) {
      Rprintf("No stream created yet!");
      return NULL;
    }
    rn = RngStream_RandU01(current_stream);
    return &rn;
  }

  void test_rstream2(double * x) {
    //unsigned long seed[6] = {12345,12345,12345,12345,12345,12345};
    //RngStream_SetPackageSeed(seed);
    Rng * s1 = new Rng("s1");
    Rng * s2 = new Rng("s2");
    // s1->set();
    x[0]=WithRNG(s1,rexp(1.0));
    x[1]=WithRNG(s2,rexp(1.0));
    s1->nextSubstream();
    x[2]=rexp(1.0);
    // current_stream.sample = NULL;
    // current_stream.state = NULL;
    delete s1;
    delete s2;
  }
  
} // extern "C"

