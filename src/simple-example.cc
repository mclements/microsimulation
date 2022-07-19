#include <microsimulation.h>

namespace {

  enum state_t {Healthy,Cancer,Death};

  enum event_t {toOtherDeath, toCancer, toCancerDeath};

  typedef std::map<std::string, std::vector<double> > Report;

  class SimplePerson : public ssim::cProcess
  {
  public:
    int id;
    state_t state;
    Report report;
    SimplePerson() : id(-1) {};
    void init();
    virtual void handleMessage(const ssim::cMessage* msg);
    void reporting(string name, double value);
  };

  /**
      Initialise a simulation run for an individual
  */
  void SimplePerson::init() {
    id++;
    state = Healthy;
    double tm = R::rweibull(8.0,85.0);
    scheduleAt(tm,toOtherDeath);
    scheduleAt(R::rweibull(3.0,90.0),toCancer);
  }

  void SimplePerson::reporting(std::string name, double value)  {
    report[name].push_back(value);
  }

  /**
      Handle receiving self-messages
  */
  void SimplePerson::handleMessage(const ssim::cMessage* msg) {

    reporting("id", double(id));
    reporting("startTime", previousEventTime);
    reporting("endtime", ssim::now());
    reporting("state", double(state));
    reporting("event", double(msg->kind));

    switch(msg->kind) {

    case toOtherDeath:
    case toCancerDeath:
      ssim::Sim::stop_process();
      break;
    
    case toCancer:
      state = Cancer;
      if (R::runif(0.0,1.0) < 0.5)
	scheduleAt(ssim::now() + R::rweibull(2.0,10.0), toCancerDeath);
      break;

    default:
      REprintf("No valid kind of event\n");
      break;
    
    } // switch

    if (id % 10000 == 0) Rcpp::checkUserInterrupt();

  } // handleMessage()

  RcppExport SEXP callSimplePerson(SEXP parms) {
    SimplePerson person;
    Rcpp::RNGScope scope;
    Rcpp::List parmsl(parms);
    int n = Rcpp::as<int>(parmsl["n"]);
    for (int i = 0; i < n; i++) {
      ssim::Sim::create_process(&person);
      ssim::Sim::run_simulation();
      ssim::Sim::clear();
    }
    return Rcpp::wrap(person.report);
  }
 
} // anonymous namespace
