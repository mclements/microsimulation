#include "microsimulation.h"

#include <boost/algorithm/cxx11/iota.hpp>

namespace {

  using namespace std;
  using namespace Rcpp;
  using namespace ssim;

  // declarations

  enum hpv_t {LR_HPV,HPV_16,HPV_18,Other_HR_HPV};
  
  enum state_t {Normal, HPV, CIN1, CIN23, LocalCancer, RegionalCancer, DistantCancer, Death};
  
  enum event_t {toHPV, toCIN1, toNormal, toCIN23, toNoCIN, toLocalCancer, toRegionalCancer,
		toDistantCancer, toUtility, toUtilityChange, toOtherDeath, toCancerDeath};

  // namespace FullState {
  //   typedef boost::tuple<short,short,short,bool,double> Type;
  //   enum Fields {state, ext_grade, dx, psa_ge_3, cohort};
  //   // string names[5] = {"state","ext_grade","dx","psa_ge_3","cohort"};
  // }
  // EventReport<FullState::Type,short,double> report;
  // typedef pair<string,double> CostKey;
  // CostReport<CostKey> costs;
  SimpleReport<double> lifeHistories;
  SimpleReport<double> outParameters;
  
  bool debug = false;

  typedef boost::tuple<hpv_t,state_t,state_t> H_key;
  typedef map<H_key,NumericInterpolate> H_t;
  H_t H;

  // typedef Table<boost::tuple<double,double,int>,double> TablePrtx; // Age, DxY, G
  // typedef Table<boost::tuple<int,double,double,int>,double> TablePradt;
  // typedef Table<pair<double,double>,double> TableBiopsyCompliance;
  // typedef Table<pair<double,double>,double> TableDDD; // as per TableBiopsyCompliance
  // typedef map<int,NumericInterpolate> H_dist_t;
  // typedef map<pair<double,int>,NumericInterpolate> H_local_t;
  // TablePrtx prtxCM, prtxRP;
  // TablePradt pradt;
  // TableBiopsyCompliance tableOpportunisticBiopsyCompliance, tableFormalBiopsyCompliance;
  // TableDDD rescreen_shape, rescreen_scale, rescreen_cure;
  // NumericInterpolate interp_prob_grade7;
  // H_dist_t H_dist;
  // H_local_t H_local;
  // set<double,greater<double> > H_local_age_set;

  Rng * rngNh, * rngOther, * rngScreen, * rngTreatment;
  Rpexp rmu0;

  NumericVector parameter;
  LogicalVector bparameter;

  class cMessageByHPV : public cMessage {
  public:
    cMessageByHPV(event_t event, hpv_t hpv) : cMessage(event), hpv(hpv) { }
    hpv_t hpv;
  };

  inline bool cMessageHPVPred(const ssim::Event* e, const hpv_t hpv) {
    const cMessageByHPV * msg = dynamic_cast<const cMessageByHPV *>(e);
    return (msg != 0 && msg->hpv == hpv); 
  }

  inline bool cMessageAnyHPVPred(const ssim::Event* e) {
    return (dynamic_cast<const cMessageByHPV *>(e) != 0); 
  }
  
  class cMessageUtilityChange : public cMessage {
  public:
    cMessageUtilityChange(double change) : cMessage(toUtilityChange), change(change) { }
    double change;
  };

  class cMessageUtility : public cMessage {
  public:
    cMessageUtility(double utility) : cMessage(toUtility), utility(utility) { }
    double utility;
  };

  // utilities
  template<class T>
  T bounds(T x, T a, T b) {
    return (x<a)?a:((x>b)?b:x);
  }
  template<class T> 
  T max(T left, T right) {
    return (left<right) ? right : left;
  }

  // forward declarations
  class Person;
  
  class HPV_infection : public cProcess
  {
  public:
    hpv_t hpv;
    state_t state;
    bool ever_infected, immunity;
    double p_immunity;
    Person * person;
    HPV_infection(hpv_t hpv = LR_HPV, double p_immunity = 1.0) : 
      hpv(hpv), p_immunity(p_immunity) {}
    void init();
    virtual void handleMessage(const cMessage* msg);
    void removeEvents() {
      Sim::remove_event(boost::bind(cMessageHPVPred,_1,hpv));
    }
  };

  class Person : public cProcess 
  {
  public:
    typedef vector<HPV_infection *> infections_t;
    int id;
    double cohort, utility;
    bool dx;
    state_t state;
    infections_t infections;
    Person(const int id = 0, const double cohort = 1950) : 
      id(id), cohort(cohort), utility(1.0) { }
    void init();
    virtual void handleMessage(const cMessage* msg);
    void link(HPV_infection * infection) { infections.push_back(infection); infection->person=this; }
    state_t max_infection_state();
    void removeHPVevents() {
      // for (infections_t::iterator it = infections.begin(); it!=infections.end(); ++it)
      // 	(*it)->removeEvents();
      Sim::remove_event(cMessageAnyHPVPred);
    }
    // void add_costs(string item);
    // void scheduleUtilityChange(double at, std::string category, bool transient = true, 
    // 			       double sign = -1.0);
  };

  state_t Person::max_infection_state() {
    int max_value = (int) Normal;
    for (vector<HPV_infection *>::iterator it = infections.begin(); it!=infections.end(); ++it)
      max_value = max<int>((int) (*it)->state, max_value);
    return (state_t) max_value;
  }

  void HPV_infection::init() {
    // initialise state
    state = Normal;
    ever_infected = false;
    immunity = false;
    // determine event times
    rngNh->set();
    scheduleAt(H[H_key(hpv,Normal,HPV)].invert(-log(R::runif(0.0,1.0))),
	       new cMessageByHPV(toHPV,hpv));
  }

  void HPV_infection::handleMessage(const cMessage* msg) {

    rngNh->set();

    if (person->id < parameter["nLifeHistories"]) { // only record up to the first n individuals
      lifeHistories.record("id",person->id);
      lifeHistories.record("time",now());
      lifeHistories.record("hpv",hpv);
      lifeHistories.record("event",msg->kind);
    }
    
    switch(msg->kind) { 
    case toNormal:
      state = Normal;
      removeEvents();
      if (!immunity)
	scheduleAt(H[H_key(hpv,Normal,HPV)].invert(-log(R::runif(0.0,1.0))),
		   new cMessageByHPV(toHPV,hpv));
      break;
    case toHPV:
      state = HPV;
      ever_infected = true;
      removeEvents();
      scheduleAt(H[H_key(hpv,HPV,Normal)].invert(-log(R::runif(0.0,1.0)), now()),
		 new cMessageByHPV(toNormal,hpv));
      scheduleAt(H[H_key(hpv,HPV,CIN1)].invert(-log(R::runif(0.0,1.0)), now()),
		 new cMessageByHPV(toCIN1,hpv));
      scheduleAt(H[H_key(hpv,HPV,CIN23)].invert(-log(R::runif(0.0,1.0)), now()),
		 new cMessageByHPV(toCIN23,hpv));
      break;
    case toCIN1:
      state = CIN1;
      removeEvents();
      scheduleAt(H[H_key(hpv,CIN1,Normal)].invert(-log(R::runif(0.0,1.0)), now()),
		 new cMessageByHPV(toNormal,hpv));
      if (!immunity)
	scheduleAt(H[H_key(hpv,CIN1,HPV)].invert(-log(R::runif(0.0,1.0)), now()),
		   new cMessageByHPV(toHPV,hpv));
      scheduleAt(H[H_key(hpv,CIN1,CIN23)].invert(-log(R::runif(0.0,1.0)), now()),
		 new cMessageByHPV(toCIN23,hpv));
      break;
    case toCIN23:
      state = CIN23;
      removeEvents();
      scheduleAt(H[H_key(hpv,CIN1,Normal)].invert(-log(R::runif(0.0,1.0)), now()),
		 new cMessageByHPV(toNormal,hpv));
      if (!immunity)
	scheduleAt(H[H_key(hpv,CIN1,HPV)].invert(-log(R::runif(0.0,1.0)), now()),
		   new cMessageByHPV(toHPV,hpv));
      scheduleAt(H[H_key(hpv,CIN1,CIN23)].invert(-log(R::runif(0.0,1.0)), now()),
		 new cMessageByHPV(toCIN23,hpv));
      break;
    case toLocalCancer:
      state = LocalCancer;
      person->removeHPVevents();
      person->scheduleAt(now(), new cMessage(toLocalCancer));
      break;
    };
  }

  /** 
      Report on costs for a given item
  */
  // void Person::add_costs(string item) {
  //   costs.add(CostKey(item,cohort),now(),cost_parameters[item]);
  // }

  /**
     Schedule a transient utility change.
     Default: sign = -1
   **/
  // void Person::scheduleUtilityChange(double at, std::string category, bool transient, double sign) {
  //   scheduleAt(at, new cMessageUtilityChange(sign*utility_estimates[category]));
  //   if (transient) {
  //     scheduleAt(at + utility_duration[category], 
  // 		 new cMessageUtilityChange(-sign*utility_estimates[category]));
  //   }
  // }

  /** 
      Initialise a simulation run for an individual
  */
  void Person::init() {
    
    // declarations
    
    // initialise state variables
    dx = false;

    // determine event times
    rngNh->set();
    
    // schedule natural history events
    double aoc = rmu0.rand(R::runif(0.0,1.0));
    scheduleAt(aoc, toOtherDeath);
    
    // schedule screening events
    // rngScreen->set();
    
    rngNh->set();
    
    // // utilities
    // // (i) set initial baseline utility
    // utility = 0.98;
    // // (ii) schedule changes in the baseline by age
    // scheduleAt(20.0, new cMessageUtility(0.97));
    // scheduleAt(40.0, new cMessageUtility(0.96));
    // scheduleAt(60.0, new cMessageUtility(0.95));
    // scheduleAt(80.0, new cMessageUtility(0.91));
    // // Mark, which reference do/should we use for this?
    
    // record some parameters using SimpleReport - too many for a tuple
    if (id<parameter["nLifeHistories"]) {
      outParameters.record("id",double(id));
      outParameters.record("aoc",aoc);
      outParameters.record("cohort",cohort);
    }
  }

/** 
    Handle self-messages received
 */
void Person::handleMessage(const cMessage* msg) {
  
  // by default, use the natural history RNG
  rngNh->set();

  // declarations
  
  // reports
  // if (bparameter["full_report"])
  //   report.add(FullState::Type(state, ext_grade, dx, psa>=3.0, cohort), msg->kind, previousEventTime, age, utility);
  if (id<parameter["nLifeHistories"]) { // only record up to the first n individuals
    lifeHistories.record("id",id);
    lifeHistories.record("time",now());
    lifeHistories.record("hpv",-1.0);
    lifeHistories.record("event",msg->kind);
  }

  // handle messages by kind

  switch(msg->kind) {

  case toCancerDeath:
    // add_costs("CancerDeath"); // cost for death, should this be zero???
    if (id<parameter["nLifeHistories"]) {
      outParameters.record("dx",dx);
      outParameters.record("age_d",now());
      outParameters.record("cancer_death",1.0);
    }
    Sim::stop_simulation();
    break;

  case toOtherDeath:
    if (id<parameter["nLifeHistories"]) {
      outParameters.record("dx",dx);
      outParameters.record("age_d",now());
      outParameters.record("cancer_death",0.0);
    }
    Sim::stop_simulation();
    break;

  case toLocalCancer:
    state = LocalCancer;
    // progression?
    break;
  
  case toRegionalCancer:
    state = RegionalCancer;
    // progression?
    break;
  
  case toDistantCancer:
    state = DistantCancer;
    // progression?
    break;


  // case toUtility:
  //   {
  //     const cMessageUtility * msgUtility;
  //     if ((msgUtility = dynamic_cast<const cMessageUtility *>(msg)) != 0) { 
  // 	utility = msgUtility->utility;
  //     } else {
  // 	REprintf("Could not cast to cMessageUtility.");
  //     }
  //   } break;

  // case toUtilityChange:
  //   {
  //     const cMessageUtilityChange * msgUtilityChange;
  //     if ((msgUtilityChange = dynamic_cast<const cMessageUtilityChange *>(msg)) != 0) { 
  // 	utility += msgUtilityChange->change;
  //     } else {
  // 	REprintf("Could not cast to cMessageUtilityChange.");
  //     }
  //   } break;

  default:
    REprintf("No valid kind of event: %i\n",msg->kind);
    break;
    
  } // switch

} // handleMessage()


RcppExport SEXP callCervical(SEXP parmsIn) {

  // declarations
  Person person;
  HPV_infection lr_hpv, hpv_16, hpv_18, other_hr_hpv;

  rngNh = new Rng();
  rngOther = new Rng();
  rngScreen = new Rng();
  rngTreatment = new Rng();
  rngNh->set();

  // read in the parameters
  List parms(parmsIn);
  List tables = parms["tables"];
  parameter = parms["parameter"];
  bparameter = parms["bparameter"]; // scalar bools
  List otherParameters = parms["otherParameters"];
  debug = as<bool>(bparameter["debug"]);
  NumericVector mu0 = as<NumericVector>(otherParameters["mu0"]);
  // cost_parameters = as<NumericVector>(otherParameters["cost_parameters"]);
  // utility_estimates = as<NumericVector>(otherParameters["utility_estimates"]);
  // utility_duration = as<NumericVector>(otherParameters["utility_duration"]);

  {
    H.clear(); 
    DataFrame df_H = as<DataFrame>(tables["H"]); // age,hpv,from_state,to_state,survival
    // extract the columns
    IntegerVector hpv = df_H["from_state"];
    IntegerVector from_state = df_H["from_state"];
    IntegerVector to_state = df_H["to_state"];
    NumericVector 
      age = df_H["age"],
      survival = df_H["survival"];
    typedef pair<double,double> dpair;
    // push to the map values and set of ages
    for (int i=0; i<age.size(); ++i) {
      // H_age_set.insert(ages[i]);
      H[H_key((hpv_t) hpv[i], (state_t) from_state[i], (state_t) to_state[i])].push_back(dpair(age[i],-log(survival[i])));
    }
    // prepare the map values for lookup
    for (H_t::iterator it = H.begin(); it != H.end(); it++) 
      it->second.prepare();
  // to use: H[H_key(hpv,from_state,to_state)].invert(-log(u),now())
  }

  NumericVector cohort = as<NumericVector>(parms["cohort"]); // at present, this is the only chuck-specific data

  // set up the parameters
  double ages0[106];
  boost::algorithm::iota(ages0, ages0+106, 0.0);
  rmu0 = Rpexp(&mu0[0], ages0, 106);
  vector<double> ages(101);
  boost::algorithm::iota(ages.begin(), ages.end(), 0.0);
  ages.push_back(1.0e+6);

  // re-set the output objects
  // report.clear();
  // costs.clear();
  outParameters.clear();
  lifeHistories.clear();

  // report.discountRate = parameter["discountRate.effectiveness"];
  // report.setPartition(ages);
  // costs.discountRate = parameter["discountRate.costs"];
  // costs.setPartition(ages);

  // main loop
  for (int i = 0; i < parameter["n"]; ++i) { 
    // define infections
    int id = i+parameter["firstID"];
    lr_hpv = HPV_infection(LR_HPV);
    hpv_16 = HPV_infection(HPV_16);
    hpv_18 = HPV_infection(HPV_18);
    other_hr_hpv = HPV_infection(Other_HR_HPV);
    // define person
    person = Person(id,cohort[i]);
    // link infections to the person
    person.link(&lr_hpv);
    person.link(&hpv_16);
    person.link(&hpv_18);
    person.link(&other_hr_hpv);
    // create processes
    Sim::create_process(&person);
    Sim::create_process(&lr_hpv);
    Sim::create_process(&hpv_16);
    Sim::create_process(&hpv_18);
    Sim::create_process(&other_hr_hpv);
    // run!!
    Sim::run_simulation();
    // tidy up
    Sim::clear();
    rngNh->nextSubstream();
    rngOther->nextSubstream();
    rngScreen->nextSubstream();
    rngTreatment->nextSubstream();
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
  }

  // tidy up
  delete rngNh;
  delete rngOther;
  delete rngScreen;
  delete rngTreatment;

  // output
  // TODO: clean up these objects in C++ (cf. R)
  return List::create(
		      // _("costs") = costs.wrap(),                // CostReport
		      // _("summary") = report.wrap(),             // EventReport 
		      _("lifeHistories") = lifeHistories.wrap(), // SimpleReport<double>
		      _("parameters") = outParameters.wrap()     // SimpleReport<double>
		      );
}

} // anonymous namespace
