
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

#include "microsimulation.h"

double rweibullHR(double shape, double scale, double hr){
  return R::rweibull(shape, scale*pow(hr,1.0/shape));
}

void remove_name(string name) {
  Sim::remove_event(bind2nd(cMessageNameEq(),name));
}

void remove_kind(short kind) {
  Sim::remove_event(bind2nd(cMessageKindEq(),kind));
}

Time now() {
  return Sim::clock();
}

Time simTime() {
  return Sim::clock();
}



static RngStream current_stream;
static double rn = 0.0;

Rng::Rng(std::string n) {
  stream = RngStream_CreateStream(n.c_str());
}

Rng::~Rng() {
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
  set(); // is this useful?
  RngStream_ResetNextSubstream(stream);
}

extern "C" {

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
