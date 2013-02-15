#include "microsimulation.h"

namespace simpleExample2 {

using namespace std;

enum state_t {Healthy,Cancer,Death};

enum event_t {toOtherDeath, toCancer, toCancerDeath};

class SimplePerson : public cProcess 
{
public:
  state_t state;
  int id;
  SimplePerson(const int i = 0) : id(i) {};
  void init();
  virtual void handleMessage(const cMessage* msg);
  static EventReport<short,short,double> report;
};

  EventReport<short,short,double> SimplePerson::report;

/** 
    Initialise a simulation run for an individual
 */
void SimplePerson::init() {
  state = Healthy;
  scheduleAt(R::rweibull(8.0,85.0),toOtherDeath);
  scheduleAt(R::rweibull(3.0,90.0),toCancer);
}

/** 
    Handle receiving self-messages
 */
void SimplePerson::handleMessage(const cMessage* msg) {

  double dwellTime, pDx;

  SimplePerson::report.add(vector<short>(1,state),msg->kind,previousEventTime,now());

  switch(msg->kind) {

  case toOtherDeath: 
  case toCancerDeath: 
    Sim::stop_simulation();
    break;
    
  case toCancer:
    state = Cancer;
    if (R::runif(0.0,1.0) < 0.5)
      scheduleAt(now() + R::rweibull(2.0,10.0), toCancerDeath);
    break;
  
  default:
    REprintf("No valid kind of event\n");
    break;
    
  } // switch

} // handleMessage()

RcppExport SEXP callSimplePerson2(SEXP parms) {
  SimplePerson person;
  Rcpp::RNGScope scope;
  Rcpp::List parmsl(parms);
  int n = Rcpp::as<int>(parmsl["n"]);

  SimplePerson::report.clear();
  vector<double> ages;
  for (double age=0.0; age<=100.0; age++) {
    ages.push_back(age);
  }
  ages.push_back(1.0e+6);
  SimplePerson::report.setPartition(ages);

  for (int i = 0; i < n; i++) {
    person = SimplePerson(i);
    Sim::create_process(&person);
    Sim::run_simulation();
    Sim::clear();
  }
  return SimplePerson::report.out();
} 

} // namespace simpleExample2
