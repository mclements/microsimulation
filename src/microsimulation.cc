#include "microsimulation.h"

double rweibullHR(double shape, double scale, double hr){
  return rweibull(shape, scale*pow(hr,1.0/shape));
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
