#include "microsimulation.h"

namespace {

  using namespace std;
  using namespace ssim;

  // declarations (NB: any changes here should also be incorporated into R code)
  
  enum stage_t {Healthy=-1,LocoRegional,Metastatic,Localised,Advanced}; // two and three stages

  enum grade_t {NoGrade=-1,GleasonLE7,GleasonGE8,GleasonLE6,Gleason7};  // two and three grades
  
  enum diagnosis_t {NotDiagnosed,ClinicalDiagnosis,ScreenDiagnosis};
  
  enum event_t {toLocalised,toAdvanced,toMetastatic,toClinicalDiagnosis,
		toCancerDeath,toOtherDeath,toScreen, 
		toBiopsy,toScreenDiagnosis,toOrganised};

  enum screen_t {noScreening, randomScreen50to70, twoYearlyScreen50to70, fourYearlyScreen50to70, screen50, screen60, screen70, screenUptake, stockholm3_goteborg, stockholm3_risk_stratified};


  typedef boost::tuple<short,short,bool,double> FullState; // stage, dx, psa_ge_3, cohort
  EventReport<FullState,short,double> report;
  map<string, vector<double> > lifeHistories; 
  map<string, vector<double> > parameters;

  //namespace par {

    /* general natural history model parameters */
  double c_age_offset = 35.0, // age offset in years for PSA calculations
    c_low_grade_slope=-0.006, // dependence of chance of low grade on age
    c_grade_num=2, // use 2 or 3 grades in natural history (assumed fixed)
    c_stage_num=2; // use 2 or 3 stages in natural history
    bool c_grade_specific_psa = true; // cf. double (is this now required?)

    /* PSA growth parameters (both all-grade and grade-specific) */
    double c_intercept_mean,
      c_intercept_sd,
      c_before_mean,
      c_before_sd,
      c_after_mean[2], // array by grade
      c_after_sd[2], // array by grade
      c_noise_sd;
      
    /* hazard rate parameters (both all-grade and grade-specific) */
    double c_onset_rate,
      c_metastasis_rate,
      c_clinical_rate_baseline,
      c_clinical_rate_distant,
      c_advanced_rate; // only if THREE_STAGE_SEQUENCE

    /* biopsy parameters */
    double c_biopsy_comp=0;  // biopsy compliance [0, 1] or -1 for PLCO rates
    double c_biopsy_sens;    // biopsy sensitivity [0, 1] or -1 for population trend

    
    // input parameters
    double screeningCompliance = 0.50;
    double studyParticipation = 0.0;
    int nLifeHistories = 10, screen = 0;

  //} // namespace par


  Rng * rngNh, * rngOther, * rngScreen;
  vector<short> stateTuple;
  Rpexp rmu0;

#define KindPred(kind) cMessageKindEq kind##Pred(kind)
#define RemoveKind(kind) Sim::remove_event(& kind##Pred)
  //cMessageKindEq toClinicalDiagnosisPred(toClinicalDiagnosis);
  KindPred(toClinicalDiagnosis);
  KindPred(toMetastatic);
  KindPred(toScreen);


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
    double age_o, y_o, age_m, age_c, age_mc;
    stage_t stage, extended_stage;
    grade_t grade, extended_grade;
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
      double yt = FhcrcPerson::ymean(t)*exp(R::rnorm(0.0, c_noise_sd));
      return yt;
    }


/** 
    Initialise a simulation run for an individual
 */
void FhcrcPerson::init() {
  
  // declarations
  double t_o, y_m, age_od;

  // change state variables
  stage = extended_stage = Healthy;
  dx = NotDiagnosed;
  everPSA = previousNegativeBiopsy = organised = false;
  rngNh->set();
  age_o = sqrt(2*R::rexp(1.0)/c_onset_rate) + c_age_offset;
  t_o = age_o - c_age_offset; // convenience
  grade = (R::runif(0.0, 1.0)>=1+c_low_grade_slope*t_o) ? GleasonLE7 : GleasonGE8; // fixed
  // TODO: extended_grade
  // TODO: treatment
  // TODO: survival
  // TODO (later): extended_stage
  beta0 = R::rnorm(c_intercept_mean,c_intercept_sd); 
  beta1 = R::rnormPos(c_before_mean,c_before_sd); 
  beta2 = R::rnormPos(c_after_mean[grade],c_after_sd[grade]); 
  y_o = ymean(t_o);
  age_m = (log((beta1+beta2)*R::rexp(1.0)/c_metastasis_rate + y_o) - beta0 + beta2*t_o) / (beta1+beta2) +
    c_age_offset;
  y_m = ymean(age_m - c_age_offset);
  age_c = (log((beta1+beta2)*R::rexp(1.0)/c_clinical_rate_baseline + yo) - beta0 + beta2*t_o) / (beta1+beta2);
  age_mc = (log((beta1+beta2)*R::rexp(1.0)/c_clinical_rate_distant + ym) - beta0 + beta2*t_o) / (beta1+beta2);
  age_od = rmu0.rand();

  // schedule natural history events
  scheduleAt(age_o,toLocalised);
  scheduleAt(age_od,toOtherDeath);

  // schedule screening events
  rngScreen->set();
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
  if (R::runif(0.0,1.0)<5.0/26.0 &&
      ((screen == stockholm3_goteborg) || (screen == stockholm3_risk_stratified)) && 
      (2013.0-cohort>=50.0 && 2013.0-cohort<70.0)) {
    scheduleAt(R::runif(2013.0,2015.0) - cohort, toOrganised);
  }

  rngNh->set();

  // record some parameters
  // faster: initialise the length of the vectors and use an index
  if (id<nLifeHistories) {
    record(parameters,"beta0",beta0);
    record(parameters,"beta1",beta1);
    record(parameters,"beta2",beta2);
    record(parameters,"age_o",age_o);
    record(parameters,"age_m",age_m);
    record(parameters,"age_c",age_c);
    record(parameters,"age_mc",age_mc);
    record(parameters,"y_o",y_o);
    record(parameters,"y_m",y_m);
    record(parameters,"stage",stage);
    record(parameters,"age_od",age_od);
    record(parameters,"cohort",cohort);
  }
}

/** 
    Handle self-messages received
 */
void FhcrcPerson::handleMessage(const cMessage* msg) {

  // declarations
  double psa = y(now()-c_age_offset); // measured PSA value (cf. ymean())
  //double year = now() + cohort;

  // record information (three states, event type, start time, end time)
  report.add(FullState(stage, dx, psa>=3.0, cohort), msg->kind, previousEventTime, now());

  if (id<nLifeHistories) { // only record up to the first n rows
    record(lifeHistories,"id", (double) id);
    record(lifeHistories,"stage", (double) stage);
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
    stage = LocoRegional;
    scheduleAt(age_c,toClinicalDiagnosis);
    scheduleAt(age_m,toMetastatic);
    break;
  
  case toMetastatic:
    stage = Metastatic;
    RemoveKind(toClinicalDiagnosis);
    scheduleAt(age_mc,toClinicalDiagnosis);
    break;
    
  case toClinicalDiagnosis:
    dx = ClinicalDiagnosis;
    RemoveKind(toMetastatic); // competing events
    RemoveKind(toScreen);
    scheduleAt(now(), toBiopsy); // for reporting (assumes three biopsies per clinical diagnosis)
    scheduleAt(now(), toBiopsy);
    scheduleAt(now(), toBiopsy);
    switch(stage) {
    case LocoRegional:
      if (R::runif(0.0,1.0) < 0.5) // 50% not cured
	scheduleAt(now() + R::rweibull(2.0,10.0), toCancerDeath);
      break;
    case Metastatic:
      if (R::runif(0.0,1.0) < 0.75) // 75% not cured
	scheduleAt(now() + R::rweibull(2.0,3.0), toCancerDeath);
      break;
    default:
      REprintf("Stage not matched\n");
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
      rngScreen->set();
      if (organised) {
	switch (screen) {
	case stockholm3_goteborg:
	  scheduleAt(now() + 2.0, toScreen);
	  break;
	case stockholm3_risk_stratified:
	  if (psa<1.0)
	    scheduleAt(now() + 8.0, toScreen);
	  else 
	    scheduleAt(now() + 4.0, toScreen);
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
    switch(stage) {
    case LocoRegional:
      if (R::runif(0.0,1.0) < 0.45) // assume slightly better stage-specific cure: 55% cf. 50%
	scheduleAt(age_c + R::rweibull(2.0,10.0), toCancerDeath);
      break;
    case Metastatic:
      if (R::runif(0.0,1.0) < 0.70) // assume slightly better stage-specific cure: 30% cf. 25%
	scheduleAt(age_mc + R::rweibull(2.0,3.0), toCancerDeath);
      break;
    default:
      REprintf("Stage not matched\n");
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
      switch(stage) {
      case Healthy:
	previousNegativeBiopsy = true;
	scheduleAt(now() + 1.0, toScreen);
	break;
      case LocoRegional:
      case Metastatic:
	scheduleAt(now(), toScreenDiagnosis);
	break;
      default:
	REprintf("Stage not matched\n");
	break;
      } break;
    } break;

  default:
    REprintf("No valid kind of event: %i\n",msg->kind);
    break;
    
  } // switch

} // handleMessage()


RcppExport SEXP callFhcrcReimplementation(SEXP parms) {

  // declarations
  FhcrcPerson person;

  rngNh = new Rng();
  rngOther = new Rng();
  rngScreen = new Rng();
  rngNh->set();

  // read in the parameters
  Rcpp::List parmsl(parms);
  int n = Rcpp::as<int>(parmsl["n"]);
  nLifeHistories = Rcpp::as<int>(parmsl["nLifeHistories"]);
  screen = Rcpp::as<int>(parmsl["screen"]);
  screeningCompliance = Rcpp::as<double>(parmsl["screeningCompliance"]);
  studyParticipation = Rcpp::as<double>(parmsl["studyParticipation"]);
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
    rngScreen->nextSubstream();
  }

  // tidy up
  delete rngNh;
  delete rngOther;
  delete rngScreen;

  // output
  return Rcpp::List::create(Rcpp::Named("summary") = report.out(),
			    Rcpp::Named("lifeHistories") = Rcpp::wrap(lifeHistories),
			    Rcpp::Named("parameters") = Rcpp::wrap(parameters)
			    );
} 

} // anonymous namespace 
