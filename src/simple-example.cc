#include "microsimulation.h"
#include "rngstream-boost.hpp"
#include <boost/random/uniform_01.hpp>

namespace {

using namespace std;
using namespace ssim;

  boost::uniform_01<> U;
  boost::rngstream gen;

  double rweibull(double u, double shape, double scale) {
    double t = scale * exp(log(-log(u))/shape);
    return t;
  }

enum state_t {Healthy,Cancer,Death};

enum event_t {toOtherDeath, toCancer, toCancerDeath,toCheck};

class SimplePerson : public cProcess 
{
public:
  state_t state;
  int id;
  SimplePerson(const int i = 0) : id(i) {};
  void init();
  virtual void handleMessage(const cMessage* msg);
};

map<string, vector<double> > report;

/** 
    Initialise a simulation run for an individual
 */
void SimplePerson::init() {
  state = Healthy;
  double tm = rweibull(U(gen),8.0,85.0); 
  scheduleAt(tm,toCheck);
  scheduleAt(tm,toOtherDeath);
  scheduleAt(rweibull(U(gen),3.0,90.0),toCancer);
}

void Reporting(string name,double value)  {
  report[name].push_back(value);
}

/** 
    Handle receiving self-messages
 */
void SimplePerson::handleMessage(const cMessage* msg) {

  Reporting("id",id);
  Reporting("startTime",previousEventTime);
  Reporting("endtime", now());
  Reporting("state", state);
  Reporting("event", msg->kind);

  switch(msg->kind) {

  case toOtherDeath: 
  case toCancerDeath: 
    SimplePerson::clear();
    break;
    
  case toCancer:
    state = Cancer;
    if (U(gen) < 0.5)
      scheduleAt(now() + rweibull(U(gen),2.0,10.0), toCancerDeath);
    break;

  case toCheck:
    break;
  
  default:
    REprintf("No valid kind of event\n");
    break;
    
  } // switch

} // handleMessage()

RcppExport SEXP vofv() {
  using namespace Rcpp;
  vector< vector< short > > in, out;
  in.push_back(vector<short>(2,1));
  in.push_back(vector<short>(2,3));
  return Rcpp::List::create(transpose(in));
}


RcppExport SEXP callSimplePerson(SEXP parms) {
  SimplePerson person;
  Rcpp::RNGScope scope;
  Rcpp::List parmsl(parms);
  int n = Rcpp::as<int>(parmsl["n"]);
  report.clear();
  Sim sim;
  for (int i = 0; i < n; i++) {
    person = SimplePerson(i);
    sim.create_process(&person);
    sim.run_simulation();
    sim.clear();
  }
  return Rcpp::wrap(report);
} 
 
} // anonymous namespace
