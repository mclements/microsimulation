#include "microsimulation.h"

#include <boost/algorithm/cxx11/iota.hpp>

namespace {

  using namespace std;
  using namespace Rcpp;
  using namespace ssim;

  // declarations

  namespace base {
    enum grade_t {Gleason_le_7,Gleason_ge_8};
  }
  namespace ext {
    enum grade_t {Gleason_le_6,Gleason_7,Gleason_ge_8};
  }
  
  enum state_t {Healthy,Localised,Metastatic}; // stage?

  enum diagnosis_t {NotDiagnosed,ClinicalDiagnosis,ScreenDiagnosis};
  
  enum event_t {toLocalised,toMetastatic,toClinicalDiagnosis,toCancerDeath,toOtherDeath,toScreen, 
		toScreenInitiatedBiopsy,toClinicalDiagnosticBiopsy,toScreenDiagnosis,toOrganised,toTreatment,toCM,toRP,toRT,toADT,toUtilityChange, toUtility };

  enum screen_t {noScreening, randomScreen50to70, twoYearlyScreen50to70, fourYearlyScreen50to70, 
		 screen50, screen60, screen70, screenUptake, stockholm3_goteborg, stockholm3_risk_stratified};

  enum treatment_t {no_treatment, CM, RP, RT};

  typedef boost::tuple<short,short,short,bool,double> FullState; // (stage, ext_grade, dx, psa_ge_3, cohort)
  typedef boost::tuple<int,short,short,bool,short,double,double,double> PSArecord; // (id, state, ext_grade, organised, dx, age, ymean, y)
  //string astates[] = {"stage", "ext_grade", "dx", "psa_ge_3", "cohort"};
  EventReport<FullState,short,double> report;
  CostReport<string> costs;
  map<string, vector<double> > lifeHistories;
  map<string, vector<double> > outParameters;
  vector<PSArecord> psarecord;

  bool debug = false;

  typedef Table<boost::tuple<double,double,int>,double> TablePrtx; // Age, DxY, G
  typedef Table<boost::tuple<int,double,double,int>,double> TablePradt;
  typedef map<int,NumericInterpolate> H_dist_t;
  typedef map<pair<double,int>,NumericInterpolate> H_local_t;
  TablePrtx prtxCM, prtxRP;
  TablePradt pradt;
  NumericInterpolate interp_prob_grade7;
  H_dist_t H_dist;
  H_local_t H_local;
  set<double,greater<double> > H_local_age_set;

  Rng * rngNh, * rngOther, * rngScreen, * rngTreatment;
  Rpexp rmu0;

  NumericVector parameter;
  NumericVector cost_parameters, utility_estimates, utility_duration;
  NumericVector mubeta2, sebeta2; // otherParameters["mubeta2"] rather than as<NumericVector>(otherParameters["mubeta2"])
  int screen, nLifeHistories;
  bool includePSArecords;

  /** 
      Utilities to record information in a map<string, vector<double> > object
  */
  void record(map<string, vector<double> > & obj, const string variable, const double value) {
    obj[variable].push_back(value);
  }
  void revise(map<string, vector<double> > & obj, const string variable, const double value) {
    //*(--obj[variable].end()) = value;
    obj[variable].pop_back();
    obj[variable].push_back(value);
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

  template<class T>
  T bounds(T x, T a, T b) {
    return (x<a)?a:((x>b)?b:x);
  }
  
  class FhcrcPerson : public cProcess 
  {
  public:
    double beta0, beta1, beta2;
    double t0, y0, tm, tc, tmc;
    state_t state;
    diagnosis_t dx;
    base::grade_t grade;
    ext::grade_t ext_grade;
    treatment_t tx;
    bool adt; 
    double txhaz;
    int id;
    double cohort, utility;
    bool everPSA, previousNegativeBiopsy, organised;
    FhcrcPerson(const int id = 0, const double cohort = 1950) : 
      id(id), cohort(cohort), utility(1.0) { };
    double ymean(double t);
    double y(double t);
    void init();
    virtual void handleMessage(const cMessage* msg);
  };

  /** 
      Calculate the (geometric) mean PSA value at a given time (** NB: time = age - 35 **)
  */
  double FhcrcPerson::ymean(double t) {
    if (t<0.0) t = 0.0; // is this the correct way to handle PSA before age 35 years?
    double yt = t<t0 ? exp(beta0+beta1*t) : exp(beta0+beta1*t+beta2*(t-t0));
    return yt;
  }
      
  /** 
      Calculate the *measured* PSA value at a given time (** NB: time = age - 35 **)
  */
  double FhcrcPerson::y(double t) {
    double yt = FhcrcPerson::ymean(t)*exp(R::rnorm(0.0, sqrt(double(parameter["tau2"]))));
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
  everPSA = previousNegativeBiopsy = organised = adt = false;
  rngNh->set();
  t0 = sqrt(2*R::rexp(1.0)/parameter["g0"]);
  grade = (R::runif(0.0, 1.0)>=1+parameter["c_low_grade_slope"]*t0) ? base::Gleason_ge_8 : base::Gleason_le_7;
  beta0 = R::rnorm(parameter["mubeta0"],parameter["sebeta0"]); 
  beta1 = R::rnormPos(parameter["mubeta1"],parameter["sebeta1"]); 
  beta2 = R::rnormPos(mubeta2[grade],sebeta2[grade]); 
  y0 = ymean(t0); // depends on: t0, beta0, beta1, beta2
  tm = (log((beta1+beta2)*R::rexp(1.0)/parameter["gm"] + y0) - beta0 + beta2*t0) / (beta1+beta2);
  ym = ymean(tm);
  tc = (log((beta1+beta2)*R::rexp(1.0)/parameter["gc"] + y0) - beta0 + beta2*t0) / (beta1+beta2);
  tmc = (log((beta1+beta2)*R::rexp(1.0)/(parameter["gc"]*parameter["thetac"]) + ym) - beta0 + beta2*t0) / (beta1+beta2);
  aoc = rmu0.rand(R::runif(0.0,1.0));
  ext_grade= (grade==base::Gleason_le_7) ? 
    (R::runif(0.0,1.0)<=interp_prob_grade7.approx(beta2) ? ext::Gleason_7 : ext::Gleason_le_6) : 
    ext::Gleason_ge_8;
  tx = no_treatment;
  txhaz = -1.0;
  // schedule natural history events
  scheduleAt(t0+35.0,toLocalised);
  scheduleAt(aoc,toOtherDeath);

  // schedule screening events
  rngScreen->set();
  if (R::runif(0.0,1.0)<parameter["screeningCompliance"]) {
    switch(screen) {
    case noScreening:
      break; // no screening
    case randomScreen50to70:
      scheduleAt(R::runif(50.0,70.0),toScreen);
      break;
    case twoYearlyScreen50to70:
      for (double start = 50.0; start<=70.0; start += 2.0) {
	scheduleAt(start, toScreen);
      }
      break;
    case fourYearlyScreen50to70: // 50,54,58,62,66,70
      for (double start = 50.0; start<=70.0; start += 4.0) {
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
  if (R::runif(0.0,1.0)<parameter["studyParticipation"] &&
      ((screen == stockholm3_goteborg) || (screen == stockholm3_risk_stratified)) && 
      (2013.0-cohort>=50.0 && 2013.0-cohort<70.0)) {
    scheduleAt(R::runif(2013.0,2015.0) - cohort, toOrganised);
  }

  rngNh->set();

  // utilities
  // (i) set initial baseline utility
  utility = 0.98;
  // (ii) schedule changes in the baseline by age
  scheduleAt(20.0, new cMessageUtility(0.97));
  scheduleAt(40.0, new cMessageUtility(0.96));
  scheduleAt(60.0, new cMessageUtility(0.95));
  scheduleAt(80.0, new cMessageUtility(0.91));
  // Mark, which reference do/should we use for this?

  // record some parameters
  // faster: initialise the length of the vectors and use an index
  if (id<nLifeHistories) {
    record(outParameters,"id",double(id)); // cf. "parameter" !!
    record(outParameters,"beta0",beta0);
    record(outParameters,"beta1",beta1);
    record(outParameters,"beta2",beta2);
    record(outParameters,"t0",t0);
    record(outParameters,"tm",tm);
    record(outParameters,"tc",tc);
    record(outParameters,"tmc",tmc);
    record(outParameters,"y0",y0);
    record(outParameters,"ym",ym);
    record(outParameters,"aoc",aoc);
    record(outParameters,"cohort",cohort);
    record(outParameters,"ext_grade",ext_grade);
    record(outParameters,"age_psa",-1.0);
    record(outParameters,"pca_death",0.0);
  }
}

/** 
    Handle self-messages received
 */
void FhcrcPerson::handleMessage(const cMessage* msg) {

  // declarations
  double psa = y(now()-35.0);
  double Z = ymean(now()-35.0);
  // double age = now();
  double year = now() + cohort;

  // record information
  report.add(FullState(state, ext_grade, dx, psa>=3.0, cohort), msg->kind, previousEventTime, now(), utility);

  if (id<nLifeHistories) { // only record up to the first n rows
    record(lifeHistories,"id", (double) id);
    record(lifeHistories,"state", (double) state);
    record(lifeHistories,"ext_grade", (double) ext_grade);
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
    costs.add("DeathCost",now(),cost_parameters["DeathCost"]); // cost for death, should this be zero???
    
    if (id<nLifeHistories) {
      record(outParameters,"age_d",now());
      revise(outParameters,"pca_death",1.0);
    }
    Sim::stop_simulation();
    scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["BiopsyUtility"]));
    scheduleAt(now() + utility_duration["BiopsyUtilityDuration"], 
	       new cMessageUtilityChange(utility_estimates["BiopsyUtility"]));
    break;

  case toOtherDeath:
    costs.add("DeathCost",now(),cost_parameters["DeathCost"]); // cost for death, should this be zero???
    
    if (id<nLifeHistories) {
      record(outParameters,"age_d",now());
    }
    Sim::stop_simulation();
    break;

  case toLocalised:
    state = Localised;
    scheduleAt(tc+35.0,toClinicalDiagnosis);
    scheduleAt(tm+35.0,toMetastatic);
    break;
  
  case toMetastatic:
    costs.add("MetastaticCancerCost",now(),cost_parameters["MetastaticCancerCost"]);
    // Scheduling utilities for palliative therapy
    scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["MetastaticCancerUtilityPart1"]));
    scheduleAt(now() + utility_duration["MetastaticCancerUtilityDurationPart1"], 
	       new cMessageUtilityChange(utility_estimates["MetastaticCancerUtilityPart1"]));
    // Scheduling utilities for terminal illness at the end of palliative therapy
    scheduleAt(now() + utility_duration["MetastaticCancerUtilityDurationPart1"], 
	       new cMessageUtilityChange(-utility_estimates["MetastaticCancerUtilityPart2"]));
    scheduleAt(now() + utility_duration["MetastaticCancerUtilityDurationPart2"], 
	       new cMessageUtilityChange(utility_estimates["MetastaticCancerUtilityPart2"]));
    state = Metastatic;
    RemoveKind(toClinicalDiagnosis);
    RemoveKind(toUtility);
    scheduleAt(tmc+35.0,toClinicalDiagnosis);
    break;
    
  case toClinicalDiagnosis:
    dx = ClinicalDiagnosis;
    RemoveKind(toMetastatic); // competing events
    RemoveKind(toScreen);
    scheduleAt(now(), toClinicalDiagnosticBiopsy); // for reporting (assumes three biopsies per clinical diagnosis)
    scheduleAt(now(), toClinicalDiagnosticBiopsy);
    scheduleAt(now(), toClinicalDiagnosticBiopsy);
    scheduleAt(now(), toTreatment);
    // switch(state) {
    // case Localised:
    //   if (R::runif(0.0,1.0) < 0.5) // 50% not cured
    // 	scheduleAt(now() + R::rweibull(2.0,10.0), toCancerDeath);
    //   break;
    // case Metastatic:
    //   if (R::runif(0.0,1.0) < 0.75) // 75% not cured
    // 	scheduleAt(now() + R::rweibull(2.0,3.0), toCancerDeath);
    //   break;
    // default:
    //   REprintf("State not matched\n");
    //   break;
    // }
    break;

  case toOrganised:
    organised = true;
    RemoveKind(toScreen); // remove other screens
    scheduleAt(now(), toScreen); // now start organised screening
    break;

  case toScreen:
    if (includePSArecords)
      psarecord.push_back(PSArecord(id, state, ext_grade, organised, dx, now(), psa, Z));
    costs.add("InvitationCost",now(),cost_parameters["InvitationCost"]);
    scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["InvitationUtility"]));
    scheduleAt(now() + utility_duration["InvitationUtilityDuration"], 
	       new cMessageUtilityChange(utility_estimates["InvitationUtility"]));
    if (organised) {
	costs.add("FormalPSACost",now(),cost_parameters["FormalPSACost"]); //Some formal screening scenarios don't seem to be set to organised
	scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["FormalPSAUtility"]));
	scheduleAt(now() + utility_duration["FormalPSAUtilityDuration"], 
		   new cMessageUtilityChange(utility_estimates["FormalPSAUtility"]));

    } else {
      costs.add("OpportunisticPSACost",now(),cost_parameters["OpportunisticPSACost"]); //cost for opportunistic PSA test
      scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["OpportunisticPSAUtility"]));
      scheduleAt(now() + utility_duration["OpportunisticPSAUtilityDuration"], 
		   new cMessageUtilityChange(utility_estimates["OpportunisticPSAUtility"]));

    }
    if (!everPSA) {
      if (id<nLifeHistories) {
	revise(outParameters,"age_psa",now());
      }
      everPSA = true;
    } 
    if (psa>=parameter["psaThreshold"]) {
      scheduleAt(now(), toScreenInitiatedBiopsy); // immediate biopsy
    } else { // re-screening schedules
      rngScreen->set();
      if (organised) {
	switch (screen) {
	case stockholm3_goteborg:
	  if (psa<1.0)
	    scheduleAt(now() + 4.0, toScreen);
	  else
	    scheduleAt(now() + 2.0, toScreen);
	  break;
	case stockholm3_risk_stratified:
	  costs.add("FormalPSABiomarkerCost",now(),cost_parameters["FormalPSABiomarkerCost"]); //cost for opportunistic PSA test
	  scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["FormalPSAUtility"]));
	  scheduleAt(now() + utility_duration["FormalPSAUtilityDuration"], 
		     new cMessageUtilityChange(utility_estimates["FormalPSAUtility"]));

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
	  double age_screen;
	  if (psa<3)
	    age_screen = now() + R::rweibull(1.16,4.79);
	  else if (psa<4)
	    age_screen = now() + R::rweibull(0.913,2.94);
	  else
	    age_screen = now() + R::rweibull(0.335,0.188);
	  // what if age_screen >= 70.0?
	  scheduleAt(age_screen, toScreen);
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
    scheduleAt(now(), toTreatment);
    // switch(state) {
    // case Localised:
    //   if (R::runif(0.0,1.0) < 0.45) // assume slightly better stage-specific cure: 55% cf. 50%
    // 	scheduleAt(tc + 35.0 + R::rweibull(2.0,10.0), toCancerDeath);
    //   break;
    // case Metastatic:
    //   if (R::runif(0.0,1.0) < 0.70) // assume slightly better stage-specific cure: 30% cf. 25%
    // 	scheduleAt(tmc + 35.0 + R::rweibull(2.0,3.0), toCancerDeath);
    //   break;
    // default:
    //   REprintf("State not matched\n");
    //   break;
    // }
    break;

  // assumes that biopsies are 100% accurate

  // record additional biopsies for clinical diagnoses
  case toClinicalDiagnosticBiopsy:
    costs.add("BiopsyCost",now(),cost_parameters["BiopsyCost"]);
    scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["BiopsyUtility"]));
    scheduleAt(now() + utility_duration["BiopsyUtilityDuration"], 
	       new cMessageUtilityChange(utility_estimates["BiopsyUtility"]));
    break;

  case toScreenInitiatedBiopsy:
    costs.add("BiopsyCost",now(),cost_parameters["BiopsyCost"]);
    scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["BiopsyUtility"]));
    scheduleAt(now() + utility_duration["BiopsyUtilityDuration"], 
	       new cMessageUtilityChange(utility_estimates["BiopsyUtility"]));

      switch(state) {
      case Healthy:
	previousNegativeBiopsy = true;
      	if (now() < 70.0) scheduleAt(now() + 1.0, toScreen);
      	break;
      case Localised:
      case Metastatic:
	scheduleAt(now(), toScreenDiagnosis);
	break;
      default:
      	REprintf("State not matched\n");
    	break;
    } break;

  case toTreatment: {
    rngTreatment->set();
    TablePrtx::key_type key = 
      TablePrtx::key_type(bounds<double>(now(),50.0,79.0),
    			  bounds<double>(year,1973.0,2004.0),
    			  int(grade));
    double pCM = prtxCM(key);
    double pRP = prtxRP(key);
    double u = R::runif(0.0,1.0);
    tx = (u<pCM)?CM:((u<pCM+pRP)?RP:RT);
    if (debug) Rprintf("Age=%3.0f, DxY=%4.0f, stage=%i, grade=%i, tx=%d, u=%8.6f, pCM=%8.6f, pRP=%8.6f\n",now(),year,state,grade,tx,u,pCM,pRP);
    if (tx == CM)
      scheduleAt(now(), toCM);
    if (tx == RP)
      scheduleAt(now(), toRP);
    if (tx == RT)
      scheduleAt(now(), toRT);
    // check for ADT
    double pADT = 
      pradt(TablePradt::key_type(tx,
				 bounds<double>(now(),50,79),
				 bounds<double>(year,1973,2004),
				 grade));
    u = R::runif(0.0,1.0);
    if (u < pADT)  {
      adt = true;
      scheduleAt(now(), toADT);
    }
    if (debug) Rprintf("adt=%d, u=%8.6f, pADT=%8.6f\n",adt,u,pADT);
    // reset the stream
    rngNh->set();
    // calculate survival
    txhaz = (state == Localised && (tx == RP || tx == RT)) ? 0.62 : 1.0;
    double sxbenefit = 1;
    u = R::runif(0.0,1.0);
    u = pow(u,1/parameter["c_baseline_specific"]); // global improvement to baseline survival
    double lead_time = (tc+35.0) - now();
    double txbenefit = exp(log(txhaz)+log(double(parameter["c_txlt_interaction"]))*lead_time);
    if (debug) Rprintf("lead_time=%f, txbenefit=%f, u=%f, ustar=%f\n",lead_time,txbenefit,u,pow(u,1/(txbenefit*sxbenefit)));
    u = pow(u,1/(txbenefit*sxbenefit));
    // TODO: calculate survival
    double age_cancer_death = -1.0;
    if (state == Localised)
      age_cancer_death = tc + 35.0 + H_local[H_local_t::key_type(*H_local_age_set.lower_bound(bounds<double>(now(),50.0,80.0)),grade)].invert(-log(u));
    if (state == Metastatic)
      age_cancer_death = tmc + 35.0 + H_dist[grade].invert(-log(u));
    scheduleAt(age_cancer_death, toCancerDeath);
    if (debug) Rprintf("SurvivalTime=%f, u=%f\n",age_cancer_death -now(), u);
  } break;

  case toRP:
    costs.add("ProstatectomyCost",now(),cost_parameters["ProstatectomyCost"]);
    // Scheduling utilities for the first 2 months after procedure
    scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["ProstatectomyUtilityPart1"]));
    scheduleAt(now() + utility_duration["ProstatectomyUtilityDurationPart1"], 
	       new cMessageUtilityChange(utility_estimates["ProstatectomyUtilityPart1"]));
    // Scheduling utilities for the first 3-12 months after procedure
    scheduleAt(now() + utility_duration["ProstatectomyUtilityDurationPart1"], 
	       new cMessageUtilityChange(-utility_estimates["ProstatectomyUtilityPart2"]));
    scheduleAt(now() + utility_duration["ProstatectomyUtilityDurationPart2"],
	       new cMessageUtilityChange(utility_estimates["ProstatectomyUtilityPart2"]));
    break;

  case toRT:
    costs.add("RadiationTherapyCost",now(),cost_parameters["RadiationTherapyCost"]);
    // Scheduling utilities for the first 2 months after procedure
    scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["RadiationTherapyUtilityPart1"]));
    scheduleAt(now() + utility_duration["RadiationTherapyUtilityDurationPart1"], 
	       new cMessageUtilityChange(utility_estimates["RadiationTherapyUtilityPart1"]));
    // Scheduling utilities for the first 3-12 months after procedure
    scheduleAt(now() + utility_duration["RadiationTherapyUtilityDurationPart1"], 
	       new cMessageUtilityChange(-utility_estimates["RadiationTherapyUtilityPart2"]));
    scheduleAt(now() + utility_duration["RadiationTherapyUtilityDurationPart2"],
	       new cMessageUtilityChange(utility_estimates["RadiationTherapyUtilityPart2"]));
    break;

  case toCM:
    costs.add("ActiveSurveillanceCost",now(),cost_parameters["ActiveSurveillanceCost"]);
    scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["ActiveSurveillanceUtility"]));
    scheduleAt(now() + utility_duration["ActiveSurveillanceUtilityDuration"], 
	       new cMessageUtilityChange(utility_estimates["ActiveSurveillanceUtility"]));

    break; 

  case toADT:
    // costs & utlities??
    break;

  case toUtility:
    {
      const cMessageUtility * msgUtility;
      if ((msgUtility = dynamic_cast<const cMessageUtility *>(msg)) != 0) { 
	utility = msgUtility->utility;
      } else {
	REprintf("Could not cast to cMessageUtility.");
      }
    } break;

  case toUtilityChange:
    {
      const cMessageUtilityChange * msgUtilityChange;
      if ((msgUtilityChange = dynamic_cast<const cMessageUtilityChange *>(msg)) != 0) { 
	utility += msgUtilityChange->change;
      } else {
	REprintf("Could not cast to cMessageUtilityChange.");
      }
    } break;

  default:
    REprintf("No valid kind of event: %i\n",msg->kind);
    break;
    
  } // switch

} // handleMessage()


RcppExport SEXP callFhcrc(SEXP parmsIn) {

  // declarations
  FhcrcPerson person;

  rngNh = new Rng();
  rngOther = new Rng();
  rngScreen = new Rng();
  rngTreatment = new Rng();
  rngNh->set();

  // read in the parameters
  List parms(parmsIn);
  List tables = parms["tables"];
  parameter = parms["parameter"];
  List otherParameters = parms["otherParameters"];
  mubeta2 = as<NumericVector>(otherParameters["mubeta2"]);
  sebeta2 = as<NumericVector>(otherParameters["sebeta2"]);
  NumericVector mu0 = as<NumericVector>(otherParameters["mu0"]);
  int n = as<int>(parms["n"]);
  includePSArecords = as<bool>(parms["includePSArecords"]);
  int firstId = as<int>(parms["firstId"]);
  interp_prob_grade7 = 
    NumericInterpolate(as<DataFrame>(tables["prob_grade7"]));
  prtxCM = TablePrtx(as<DataFrame>(tables["prtx"]),
			       "Age","DxY","G","CM"); // NB: Grade is now {0,1} coded cf {1,2}
  prtxRP = TablePrtx(as<DataFrame>(tables["prtx"]),
			       "Age","DxY","G","RP");
  pradt = TablePradt(as<DataFrame>(tables["pradt"]),"Tx","Age","DxY","Grade","ADT");

  H_dist.clear();
  DataFrame df_survival_dist = as<DataFrame>(tables["survival_dist"]); // Grade,Time,Survival
  DataFrame df_survival_local = as<DataFrame>(tables["survival_local"]); // Age,Grade,Time,Survival
  // extract the columns from the survival_dist data-frame
  IntegerVector sd_grades = df_survival_dist["Grade"];
  NumericVector 
    sd_times = df_survival_dist["Time"],
    sd_survivals = df_survival_dist["Survival"];
  typedef pair<double,double> dpair;
  for (int i=0; i<sd_grades.size(); ++i) 
    H_dist[sd_grades[i]].push_back(dpair(sd_times[i],-log(sd_survivals[i])));
  for (H_dist_t::iterator it_sd = H_dist.begin(); it_sd != H_dist.end(); it_sd++) 
    it_sd->second.prepare();
  // now we can use: H_dist[grade].invert(-log(u))
  H_local.clear();
  // extract the columns from the data-frame
  IntegerVector sl_grades = df_survival_local["Grade"];
  NumericVector 
    sl_ages = df_survival_local["Age"],
    sl_times = df_survival_local["Time"],
    sl_survivals = df_survival_local["Survival"];
  // push to the map values and set of ages
  for (int i=0; i<sl_grades.size(); ++i) {
    H_local_age_set.insert(sl_ages[i]);
    H_local[H_local_t::key_type(sl_ages[i],sl_grades[i])].push_back
      (dpair(sl_times[i],-log(sl_survivals[i])));
  }
  // prepare the map values for lookup
  for (H_local_t::iterator it_sl = H_local.begin(); 
       it_sl != H_local.end(); 
       it_sl++) 
    it_sl->second.prepare();
  // now we can use: H_local[H_local_t::key_type(*H_local_age_set.lower_bound(age),grade)].invert(-log(u))

  if (debug) {
    Rprintf("SurvTime: %f\n",exp(-H_local[H_local_t::key_type(65.0,0)].approx(63.934032)));
    Rprintf("SurvTime: %f\n",H_local[H_local_t::key_type(*H_local_age_set.lower_bound(65.0),0)].invert(-log(0.5)));
    Rprintf("SurvTime: %f\n",exp(-H_dist[0].approx(5.140980)));
    Rprintf("SurvTime: %f\n",H_dist[0].invert(-log(0.5)));
  }

  nLifeHistories = as<int>(otherParameters["nLifeHistories"]);
  screen = as<int>(otherParameters["screen"]);
  NumericVector cohort = as<NumericVector>(parms["cohort"]); // at present, this is the only chuck-specific data

  cost_parameters = as<NumericVector>(parms["cost_parameters"]);
  utility_estimates = as<NumericVector>(parms["utility_estimates"]);
  utility_duration = as<NumericVector>(parms["utility_duration"]);

  // set up the parameters
  double ages0[106];
  boost::algorithm::iota(ages0, ages0+106, 0.0);
  rmu0 = Rpexp(&mu0[0], ages0, 106);
  vector<double> ages(101);
  boost::algorithm::iota(ages.begin(), ages.end(), 0.0);
  ages.push_back(1.0e+6);

  // re-set the output objects
  report.clear();
  costs.clear();
  outParameters.clear();
  lifeHistories.clear();
  psarecord.clear();

  report.setPartition(ages);
  costs.setPartition(ages);

  // main loop
  for (int i = 0; i < n; ++i) {
    person = FhcrcPerson(i+firstId,cohort[i]);
    Sim::create_process(&person);
    Sim::run_simulation();
    Sim::clear();
    rngNh->nextSubstream();
    rngOther->nextSubstream();
    rngScreen->nextSubstream();
    rngTreatment->nextSubstream();
  }

  // tidy up
  delete rngNh;
  delete rngOther;
  delete rngScreen;
  delete rngTreatment;

  // output
  // TODO: clean up these objects in C++ (cf. R)
  return List::create(_("costs") = costs.out(),
		      _("summary") = report.out(),
		      _("lifeHistories") = wrap(lifeHistories),
		      _("parameters") = wrap(outParameters),
		      _("psarecord")=wrap(psarecord)
		      );
}

} // anonymous namespace
