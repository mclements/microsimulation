#include <microsimulation.h>
#include <numeric>

namespace illnessDeath {

  using namespace std;
  using namespace ssim;

  enum state_t {Healthy,Cancer};
  
  enum event_t {toOtherDeath, toCancer, toCancerDeath};

  EventReport<short,short,double> report;
  double cure, zsd; // parameters - could be static class variables
  

  double b_weibull(double mean, double a, double rr = 1.0) {
    return mean/R::gammafn(1.0+1.0/a)*pow(rr,-1.0/a);
  }

  class SimplePerson : public cProcess 
  {
  public:
    state_t state;
    int id;
    double z;
    SimplePerson(const int i = 0) : id(i) {};
    void init();
    virtual void handleMessage(const cMessage* msg);
  };
  
  /** 
      Initialise a simulation run for an individual
  */
  void SimplePerson::init() {
    state = Healthy;
    z = exp(R::rnorm(0.0,zsd));
    scheduleAt(R::rweibull(4.0,b_weibull(80.0,4.0)), toOtherDeath);
    if (R::runif(0.0,1.0)>cure)
      scheduleAt(R::rweibull(3.0,b_weibull(80.0,3.0,z)), toCancer);
  }
  
  /** 
      Handle receiving self-messages
  */
  void SimplePerson::handleMessage(const cMessage* msg) {
    
    report.add(state, msg->kind, previousEventTime, now());
    
    switch(msg->kind) {
      
    case toOtherDeath: 
    case toCancerDeath: 
      // reporting already completed: stop the simulation
      Sim::stop_simulation();
      break;
      
    case toCancer:
      state = Cancer;
      RemoveKind(toOtherDeath);
      if (R::runif(0.0,1.0) < 0.5) // cure fraction
	scheduleAt(now() + R::rweibull(1.0,10.0), toCancerDeath);
      break;
      
    default:
      REprintf("No valid kind of event\n");
      break;
      
    } // switch
    
  } // handleMessage()
  
  RcppExport SEXP callIllnessDeath(SEXP parms) {
    SimplePerson person;
    Rcpp::RNGScope scope;
    Rcpp::List parmsl(parms);
    int n = Rcpp::as<int>(parmsl["n"]);
    cure = Rcpp::as<double>(parmsl["cure"]);
    zsd = Rcpp::as<double>(parmsl["zsd"]);
    
    vector<double> ages(101);
    std::iota(ages.begin(), ages.end(), 0.0);
    ages.push_back(1.0e+6);
    report.clear();
    report.setPartition(ages);
    
    for (int i = 0; i < n; i++) {
      person = SimplePerson(i);
      Sim::create_process(&person);
      Sim::run_simulation();
      Sim::clear();
    }
    return report.wrap();
  } 
  
} // namespace illnessDeath
