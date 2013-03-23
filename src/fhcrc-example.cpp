#include "microsimulation.h"

namespace fhcrc {

  using namespace std;

  // declarations
  
  enum state_t {Healthy,Localised,Metastatic};

  enum diagnosis_t {NotDiagnosed,ClinicalDiagnosis,ScreenDiagnosis};
  
  enum event_t {toLocalised,toMetastatic,toClinicalDiagnosis,toCancerDeath,toOtherDeath,toScreen, 
		toBiopsy,toScreenDiagnosis,toOrganised};

  enum screen_t {noScreening, randomScreen50to70, twoYearlyScreen50to70, fourYearlyScreen50to70, screen50, screen60, screen70, screenUptake, stockholm3_goteborg, stockholm3_risk_stratified};

  double tau2 = 0.0829,
    g0=0.0005, gm=0.0004, gc=0.0015, 
    thetac=19.1334,
    mubeta0=-1.6094, varbeta0=0.0568,
    mubeta1=0.02, varbeta1=0.0019,
    mubeta2=0.1094,varbeta2=0.0237;
  double screeningCompliance = 0.50;
  int nLifeHistories = 10, screen = 0;
  Rng * rngNh, * rngOther;
  vector<short> stateTuple;
  Rpexp rmu0;

  //#define KindPred(kind) cMessageKindEq kind##Pred = cMessageKindEq(kind)
#define RemoveKind(KIND) remove_kind(KIND)
  // #define RemoveKind(kind) Sim::remove_event(& kind##Pred)
  // KindPred(toClinicalDiagnosis);
  // KindPred(toMetastatic);
  // KindPred(toScreen);

// #define KindPred2(KIND) class KIND##Pred : public ssim::EventPredicate { \
//  public: \
//   bool operator()(const ssim::Event* e)  { \
//     const cMessage * msg = dynamic_cast<const cMessage *>(e); \
//     return (msg != 0 && msg->kind == KIND);\
//   }; \
// };


  EventReport<short,short,double> report;
  map<string, vector<double> > lifeHistories; 
  map<string, vector<double> > parameters;

  /** 
      Utility to record information in a map<string vector<double> > object
  */
  void record(map<string, vector<double> > & obj, string variable, double value) {
    obj[variable].push_back(value);
  }
  
  // all cause mortality rates by single year of age from age 0 could be a parameter input
  double mu0[] = {0.00219, 0.000304, 5.2e-05, 0.000139, 0.000141, 3.6e-05, 7.3e-05, 
		  0.000129, 3.8e-05, 0.000137, 6e-05, 8.1e-05, 6.1e-05, 0.00012, 
		  0.000117, 0.000183, 0.000185, 0.000397, 0.000394, 0.000585, 0.000448, 
		  0.000696, 0.000611, 0.000708, 0.000659, 0.000643, 0.000654, 0.000651, 
		  0.000687, 0.000637, 0.00063, 0.000892, 0.000543, 0.00058, 0.00077, 
		  0.000702, 0.000768, 0.000664, 0.000787, 0.00081, 0.000991, 9e-04, 
		  0.000933, 0.001229, 0.001633, 0.001396, 0.001673, 0.001926, 0.002217, 
		  0.002562, 0.002648, 0.002949, 0.002729, 0.003415, 0.003694, 0.004491, 
		  0.00506, 0.004568, 0.006163, 0.006988, 0.006744, 0.00765, 0.007914, 
		  0.009153, 0.010231, 0.011971, 0.013092, 0.013839, 0.015995, 0.017693, 
		  0.018548, 0.020708, 0.022404, 0.02572, 0.028039, 0.031564, 0.038182, 
		  0.042057, 0.047361, 0.05315, 0.058238, 0.062619, 0.074934, 0.089776, 
		  0.099887, 0.112347, 0.125351, 0.143077, 0.153189, 0.179702, 0.198436, 
		  0.240339, 0.256215, 0.275103, 0.314157, 0.345252, 0.359275, 0.41768, 
		  0.430279, 0.463636, 0.491275, 0.549738, 0.354545, 0.553846, 0.461538, 
		  0.782609};

  class FhcrcPerson : public cProcess 
  {
  public:
    double beta0, beta1, beta2;
    double t0, y0, tm, tc, tmc;
    state_t state;
    diagnosis_t dx;
    int id;
    double cohort;
    bool everPSA, previousNegativeBiopsy, organised;
    FhcrcPerson(const int i = 0, const double coh = 1950) : id(i), cohort(coh) { };
    double ymean(double t);
    double y(double t);
    void init();
    virtual void handleMessage(const cMessage* msg);
  };

  

  /** 
      Calculate the (geometric) mean PSA value at a given time (** NB: time = age - 35 **)
  */
  double FhcrcPerson::ymean(double t) {
    if (t<0) t = 0; // is this the correct way to handle PSA before age 35 years?
    double yt = t<t0 ? exp(beta0+beta1*t) : exp(beta0+beta1*t+beta2*(t-t0));
    return yt;
  }
      
  /** 
      Calculate the *measured* PSA value at a given time (** NB: time = age - 35 **)
  */
  double FhcrcPerson::y(double t) {
      double yt = FhcrcPerson::ymean(t)*exp(R::rnorm(0.0, sqrt(tau2)));
      return yt;
    }

/** 
    Initialise a simulation run for an individual
 */
void FhcrcPerson::init() {
  
  // declarations
  double ym, aoc;

  // change state variables
  state = Healthy;
  dx = NotDiagnosed;
  everPSA = previousNegativeBiopsy = organised = false;
  rngNh->set();
  beta0 = R::rnorm(mubeta0,sqrt(varbeta0)); 
  beta1 = R::rnormPos(mubeta1,sqrt(varbeta1)); 
  beta2 = R::rnormPos(mubeta2,sqrt(varbeta2)); 
  t0 = sqrt(2*R::rexp(1.0)/g0);
  y0 = ymean(t0);
  tm = (log((beta1+beta2)*R::rexp(1.0)/gm + y0) - beta0 + beta2*t0) / (beta1+beta2);
  ym = ymean(tm);
  tc = (log((beta1+beta2)*R::rexp(1.0)/gc + y0) - beta0 + beta2*t0) / (beta1+beta2);
  tmc = (log((beta1+beta2)*R::rexp(1.0)/(gc*thetac) + ym) - beta0 + beta2*t0) / (beta1+beta2);
  aoc = rmu0.rand();

  // schedule natural history events
  scheduleAt(t0+35.0,toLocalised);
  scheduleAt(aoc,toOtherDeath);

  // schedule screening events
  rngOther->set();
  if (R::runif(0.0,1.0)<screeningCompliance) {
    switch(screen) {
    case noScreening:
      break; // no screening
    case randomScreen50to70:
      scheduleAt(R::runif(50.0,70.0),toScreen);
      break;
    case twoYearlyScreen50to70:
      for (double start = 50.0; start<=70.0; start = start + 2.0) {
	scheduleAt(start, toScreen);
      }
      break;
    case fourYearlyScreen50to70:
      for (double start = 50.0; start<=70.0; start = start + 4.0) {
	scheduleAt(start, toScreen);
      }
      break;
    case screen50:
      scheduleAt(50.0,toScreen);
      break;
    case screen60:
      scheduleAt(60.0,toScreen);
      break;
    case screen70:
      scheduleAt(70.0,toScreen);
      break;
    case stockholm3_goteborg:
    case stockholm3_risk_stratified:
    case screenUptake: {
      // screening participation increases with time
      if (1995.0 - cohort < 50.0) 
	scheduleAt(50.0 + R::rweibull(2.0,10.0), toScreen);
      else
	scheduleAt(1995.0 - cohort + R::rweibull(2.0,10.0), toScreen);
      // for re-screening patterns: see handleMessage()
    } break;
    default:
      REprintf("Screening not matched: %i\n",screen);
      break;
    }
  }
  if (R::runif(0.0,1.0)<9.0/26.0 &&
      ((screen == stockholm3_goteborg) || (screen == stockholm3_risk_stratified)) && 
      (2013.0-cohort>=50.0 && 2013.0-cohort<70.0)) {
    scheduleAt(R::runif(2013.0,2015.0) - cohort, toOrganised);
  }

  // record some parameters
  // faster: initialise the length of the vectors and use an index
  if (id<nLifeHistories) {
    record(parameters,"beta0",beta0);
    record(parameters,"beta1",beta1);
    record(parameters,"beta2",beta2);
    record(parameters,"t0",t0);
    record(parameters,"tm",tm);
    record(parameters,"tc",tc);
    record(parameters,"tmc",tmc);
    record(parameters,"y0",y0);
    record(parameters,"ym",ym);
    record(parameters,"aoc",aoc);
    record(parameters,"cohort",cohort);
  }
}

/** 
    Handle self-messages received
 */
void FhcrcPerson::handleMessage(const cMessage* msg) {

  // declarations
  double psa = y(now()-35.0);
  //double year = now() + cohort;

  // record information (three states, event type, start time, end time)
  report.add(state, dx, psa>=3.0 ? 1 : 0, short(cohort), msg->kind, previousEventTime, now());

  if (id<nLifeHistories) { // only record up to the first n rows
    record(lifeHistories,"id", (double) id);
    record(lifeHistories,"state", (double) state);
    record(lifeHistories,"dx", (double) dx);
    record(lifeHistories,"event", (double) msg->kind);
    record(lifeHistories,"begin", previousEventTime);
    record(lifeHistories,"end", now());
    record(lifeHistories,"psa", psa);
  }

  // by default, use the natural history RNG
  rngNh->set();

  // handle messages by kind

  switch(msg->kind) {

  case toCancerDeath: 
  case toOtherDeath: 
    Sim::stop_simulation();
    break;

  case toLocalised:
    state = Localised;
    scheduleAt(tc+35.0,toClinicalDiagnosis);
    scheduleAt(tm+35.0,toMetastatic);
    break;
  
  case toMetastatic:
    state = Metastatic;
    RemoveKind(toClinicalDiagnosis);
    scheduleAt(tmc+35.0,toClinicalDiagnosis);
    break;
    
  case toClinicalDiagnosis:
    dx = ClinicalDiagnosis;
    RemoveKind(toMetastatic); // competing events
    RemoveKind(toScreen);
    scheduleAt(now(), toBiopsy); // for reporting (assumes three biopsies per clinical diagnosis)
    scheduleAt(now(), toBiopsy);
    scheduleAt(now(), toBiopsy);
    switch(state) {
    case Localised:
      if (R::runif(0.0,1.0) < 0.5) // 50% not cured
	scheduleAt(now() + R::rweibull(2.0,10.0), toCancerDeath);
      break;
    case Metastatic:
      if (R::runif(0.0,1.0) < 0.75) // 75% not cured
	scheduleAt(now() + R::rweibull(2.0,3.0), toCancerDeath);
      break;
    default:
      REprintf("State not matched\n");
      break;
    }
    break;

  case toOrganised:
    organised = true;
    RemoveKind(toScreen); // remove other screens
    scheduleAt(now(), toScreen); // now start organised screening
    break;

  case toScreen: 
    everPSA = true;
    if (psa>=3.0) {
      scheduleAt(now(), toBiopsy); // immediate biopsy
    } else { // re-screening schedules
      rngOther->set();
      if (organised) {
	switch (screen) {
	case stockholm3_goteborg:
	  scheduleAt(now() + 2.0, toScreen);
	  break;
	case stockholm3_risk_stratified:
	  if (psa<1.0)
	    scheduleAt(now() + 6.0, toScreen);
	  else 
	    scheduleAt(now() + 2.0, toScreen);
	  break;
	default:
	  REprintf("Organised screening state not matched: %s\n",screen);
	  break;
	}
      } else  // opportunistic screening
	switch(screen) {
	case screenUptake:
	case stockholm3_goteborg:
	case stockholm3_risk_stratified:
	  double tm;
	  if (psa<3)
	    tm = now() + R::rweibull(1.16,4.79);
	  else if (psa<4)
	    tm = now() + R::rweibull(0.913,2.94);
	  else
	    tm = now() + R::rweibull(0.335,0.188);
	  scheduleAt(tm, toScreen);
	  break;
	default:
	  // do not schedule any other screens
	  break;
	}
      rngNh->set();
    }
    break;
    
  case toScreenDiagnosis:
    dx = ScreenDiagnosis;
    RemoveKind(toMetastatic); // competing events
    RemoveKind(toClinicalDiagnosis);
    RemoveKind(toScreen);
    switch(state) {
    case Localised:
      if (R::runif(0.0,1.0) < 0.45) // assume slightly better stage-specific cure: 55% cf. 50%
	scheduleAt(tc + 35.0 + R::rweibull(2.0,10.0), toCancerDeath);
      break;
    case Metastatic:
      if (R::runif(0.0,1.0) < 0.70) // assume slightly better stage-specific cure: 30% cf. 25%
	scheduleAt(tmc + 35.0 + R::rweibull(2.0,3.0), toCancerDeath);
      break;
    default:
      REprintf("State not matched\n");
      break;
    }
    break;

    // assumes that biopsies are 100% accurate
  case toBiopsy:
    switch(dx) {
    case(ScreenDiagnosis):
    case(ClinicalDiagnosis):
      // already diagnosed
      break;
    case(NotDiagnosed): 
      switch(state) {
      case Healthy:
	previousNegativeBiopsy = true;
	scheduleAt(now() + 1.0, toScreen);
	break;
      case Localised:
      case Metastatic:
	scheduleAt(now(), toScreenDiagnosis);
	break;
      default:
	REprintf("State not matched\n");
	break;
      } break;
    } break;

  default:
    REprintf("No valid kind of event: %i\n",msg->kind);
    break;
    
  } // switch

} // handleMessage()

RcppExport SEXP testrexp() {
  return Rcpp::wrap(R::rexp(100.0));
}

RcppExport SEXP callFhcrc(SEXP parms) {

  // declarations
  FhcrcPerson person;

  rngNh = new Rng();
  rngOther = new Rng();
  rngNh->set();

  // read in the parameters
  Rcpp::List parmsl(parms);
  int n = Rcpp::as<int>(parmsl["n"]);
  nLifeHistories = Rcpp::as<int>(parmsl["nLifeHistories"]);
  screen = Rcpp::as<int>(parmsl["screen"]);
  screeningCompliance = Rcpp::as<double>(parmsl["screeningCompliance"]);
  Rcpp::NumericVector cohort = Rcpp::as<Rcpp::NumericVector>(parmsl["cohort"]);

  // set up the parameters
  double ages0[106];
  iota(ages0, ages0+106, 0.0);
  rmu0 = Rpexp(mu0, ages0, 106);
  vector<double> ages(101);
  iota(ages.begin(), ages.end(), 0.0);
  ages.push_back(1.0e+6);

  // re-set the output objects
  report.clear();
  parameters.clear();
  lifeHistories.clear();

  report.setPartition(ages);

  // main loop
  for (int i = 0; i < n; i++) {
    person = FhcrcPerson(i,cohort[i]);
    Sim::create_process(&person);
    Sim::run_simulation();
    Sim::clear();
    rngNh->nextSubstream();
    rngOther->nextSubstream();
  }

  // tidy up
  delete rngNh;
  delete rngOther;

  // output
  return Rcpp::List::create(Rcpp::Named("summary") = report.out(),
			    Rcpp::Named("lifeHistories") = Rcpp::wrap(lifeHistories),
			    Rcpp::Named("parameters") = Rcpp::wrap(parameters)
			    );
} 

} // namespace fhcrc
