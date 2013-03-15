
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

#include "microsimulation.h"

double rweibullHR(double shape, double scale, double hr){
  return R::rweibull(shape, scale*pow(hr,1.0/shape));
}

void remove_name(string name) {
  Sim::remove_event(cMessageNameEq(name));
}

void remove_kind(short kind) {
  Sim::remove_event(cMessageKindEq(kind));
}

// void remove_name(string name) {
//   Sim::remove_event(bind2nd(cMessageNameEq(),name));
// }

// void remove_kind(short kind) {
//   Sim::remove_event(bind2nd(cMessageKindEq(),kind));
// }

Time now() {
  return Sim::clock();
}

Time simTime() {
  return Sim::clock();
}



static RngStream current_stream = NULL;
static RngStream default_stream;
static double rn = 0.0;

Rng::Rng(std::string s) {
  stream = RngStream_CreateStream(s.c_str());
}

Rng::~Rng() {
  if (current_stream == stream)
    current_stream = default_stream;
  RngStream_DeleteStream(&stream);
}

void Rng::set() {
  if (stream == NULL) {
    REprintf("Rng::set(): Invalid NULL pointer");
  } else {
    current_stream = stream;
  }
  return;
}

void Rng::nextSubstream() {
  //set(); // is this useful?
  RngStream_ResetNextSubstream(stream);
}

void Rng::nextSubStream() {
  nextSubstream();
}

// void Rng::unset() {
//   current_stream = default_stream;
// }


extern "C" {

void r_create_current_stream()
{
  default_stream = RngStream_CreateStream("default stream");
  current_stream = default_stream;
  // current_stream = (RngStream) malloc (sizeof (struct RngStream_InfoState));
  
  // if (current_stream == NULL) {
  //   error("r_create_current_stream: No more memory\n\n");
  // }
  // current_stream->name=(char *)"";
  return;
}

void r_remove_current_stream()
{
  // free(current_stream);
  RngStream_DeleteStream(&default_stream);
  return;
}

  void r_set_user_random_seed(int * inseed) {
    unsigned long seed[6];
    for(int i=0; i<6; i++) {
      seed[i] = (unsigned long)inseed[i];
    }
    RngStream_SetSeed (default_stream, seed);
  }

  void r_get_user_random_seed(int * outseed) {
    unsigned long seed[6];
    RngStream_GetState (default_stream, seed);
    for(int i=0; i<6; i++) {
      outseed[i] = (int)seed[i];
    }
  }


double *user_unif_rand ()
{
  if (!current_stream) {
    REprintf("user_unif_rand(): No stream created yet!");
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
  x[0]=WithRNG(s1,R::rexp(1.0));
  x[1]=WithRNG(s2,R::rexp(1.0));
  s1->nextSubstream();
  x[2]=R::rexp(1.0);
  delete s1;
  delete s2;
}

} // extern "C"


double discountedInterval(double start, double end, double discountRate) {
  if (discountRate == 0.0) return end - start;
  //else if (start == 0.0) return (1.0 - (1.0+discountRate)^(-end)) / log(1.0+discountRate);
  else return (pow(1.0+discountRate,-start) - pow(1.0+discountRate,-end)) / log(1.0+discountRate);
}

namespace R {
  double rnormPos(double mean, double sd) {
    double x;
    while ((x=R::rnorm(mean,sd))<0.0) { }
    return x;
  }
}
