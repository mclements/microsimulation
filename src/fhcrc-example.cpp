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
  
  enum event_t {toLocalised,toMetastatic,toClinicalDiagnosis,toCancerDeath,toOtherDeath,toScreen,toBiopsyFollowUpScreen, 
		toScreenInitiatedBiopsy,toClinicalDiagnosticBiopsy,toScreenDiagnosis,toOrganised,toTreatment,toCM,toRP,toRT,toADT,toUtilityChange, toUtility };

  enum screen_t {noScreening, randomScreen50to70, twoYearlyScreen50to70, fourYearlyScreen50to70, 
		 screen50, screen60, screen70, screenUptake, stockholm3_goteborg, stockholm3_risk_stratified};

  enum treatment_t {no_treatment, CM, RP, RT};

  namespace FullState {
    typedef boost::tuple<short,short,short,bool,double> Type;
    enum Fields {state, ext_grade, dx, psa_ge_3, cohort};
    // string names[5] = {"state","ext_grade","dx","psa_ge_3","cohort"};
  }
  // namespace PSArecord {
  //   typedef boost::tuple<int,short,short,bool,short,double,double,double> Type;
  //   enum Fields {id,state,ext_grade,organised,dx,age,ymean,beta0,beta1,beta2,t0};
  // }
  namespace LifeHistory {
    typedef boost::tuple<int,short,short,int,short,double,double,double,double> Type;
    enum Fields {id,state,ext_grade,dx,event,begin,end,year,psa};
  }

  template<class T = double>
  class SimpleReport {
  public:
    typedef map<string,vector<T> > Map;
    void record(string field, T value) {
      _data[field].push_back(value);
    }
    void revise(string field, T value) {
      _data[field].pop_back();
      _data[field].push_back(value);
    }
    void clear() { _data.clear(); }
    SEXP wrap() {
      return Rcpp::wrap(_data);
    }
    void append(SimpleReport<T> & obj) {
      for(typename Map::iterator it = obj._data.begin(); it != obj._data.end(); ++it) {
	_data[it->first].insert(_data[it->first].end(), it->second.begin(), it->second.end());
      }
    }
    Map _data;
  };

  RcppExport SEXP rllogis_(SEXP shape, SEXP scale) {
    RNGScope scope;
    return wrap(R::rllogis(as<double>(shape),as<double>(scale)));
  }

  //typedef boost::tuple<short,short,short,bool,double> FullState; // (stage, ext_grade, dx, psa_ge_3, cohort)
  //typedef boost::tuple<int,short,short,bool,short,double,double,double> PSArecord; // (id, state, ext_grade, organised, dx, age, ymean, y)
  EventReport<FullState::Type,short,double> report;
  CostReport<string> costs;
  vector<LifeHistory::Type> lifeHistories;
  SimpleReport<double> outParameters;
  SimpleReport<double> psarecord;

  bool debug = false;

  typedef Table<boost::tuple<double,double,int>,double> TablePrtx; // Age, DxY, G
  typedef Table<boost::tuple<int,double,double,int>,double> TablePradt;
  typedef Table<pair<double,double>,double> TableBiopsyCompliance;
  typedef Table<pair<double,double>,double> TableDDD; // as per TableBiopsyCompliance
  typedef map<int,NumericInterpolate> H_dist_t;
  typedef map<pair<double,int>,NumericInterpolate> H_local_t;
  TablePrtx prtxCM, prtxRP;
  TablePradt pradt;
  TableBiopsyCompliance tableBiopsyCompliance;
  TableDDD rescreen_shape, rescreen_scale, rescreen_cure;
  NumericInterpolate interp_prob_grade7;
  H_dist_t H_dist;
  H_local_t H_local;
  set<double,greater<double> > H_local_age_set;

  Rng * rngNh, * rngOther, * rngScreen, * rngTreatment;
  Rpexp rmu0;

  NumericVector parameter;

  // read in the parameters
  NumericVector cost_parameters, utility_estimates, utility_duration;
  NumericVector mubeta2, sebeta2; // otherParameters["mubeta2"] rather than as<NumericVector>(otherParameters["mubeta2"])
  int screen, nLifeHistories;
  bool includePSArecords;

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
      double pscreening = cohort>=1932.0 ? 0.9 : 0.9-(1932.0 - cohort)*0.01; 
      double shape = 3.81;
      double scale = cohort>1960.0 ? 19.4 : 19.4*17.0/(1977.0 - cohort);
      double start_age = cohort>1960.0 ? 35.0 : 1995.0 - cohort;
      double first_screen = start_age+R::rllogis(shape,scale);
      double u = R::runif(0.0,1.0);
      if (u<pscreening)
	scheduleAt(first_screen, toScreen);
      // Rprintf("(cohort=%f,pscreening=%f,u=%f,start_age=%f,first_screen=%f)\n",cohort,pscreening,u,start_age,first_screen);
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

  // record some parameters using SimpleReport - too many for a tuple
  if (id<nLifeHistories) {
    outParameters.record("id",double(id));
    outParameters.record("beta0",beta0);
    outParameters.record("beta1",beta1);
    outParameters.record("beta2",beta2);
    outParameters.record("t0",t0);
    outParameters.record("tm",tm);
    outParameters.record("tc",tc);
    outParameters.record("tmc",tmc);
    outParameters.record("y0",y0);
    outParameters.record("ym",ym);
    outParameters.record("aoc",aoc);
    outParameters.record("cohort",cohort);
    outParameters.record("ext_grade",ext_grade);
    outParameters.record("age_psa",-1.0);
    outParameters.record("pca_death",0.0);
  }
}

/** 
    Handle self-messages received
 */
void FhcrcPerson::handleMessage(const cMessage* msg) {
  
  // declarations
  double psa = y(now()-35.0);
  double Z = ymean(now()-35.0);
  double age = now();
  double year = age + cohort;
  double compliance;

  // record information
  report.add(FullState::Type(state, ext_grade, dx, psa>=3.0, cohort), msg->kind, previousEventTime, age, utility);

  if (id<nLifeHistories) { // only record up to the first n individuals
    lifeHistories.push_back(LifeHistory::Type(id,state,ext_grade,dx,msg->kind,previousEventTime,age,year,psa));
  }
  // by default, use the natural history RNG
  rngNh->set();

  // handle messages by kind

  switch(msg->kind) {

  case toCancerDeath:
    costs.add("DeathCost",now(),cost_parameters["DeathCost"]); // cost for death, should this be zero???
    // change in utility for terminal illness at the end of palliative therapy
    // utility -= utility_duration["PalliativeDuration"]*(utility_estimates["MetastaticCancerUtility"]-utility_estimates["PalliativeUtility"]);
    if (id<nLifeHistories) {
      outParameters.record("age_d",now());
      outParameters.revise("pca_death",1.0);
    }
    Sim::stop_simulation();
    break;

  case toOtherDeath:
    costs.add("DeathCost",now(),cost_parameters["DeathCost"]); // cost for death, should this be zero???
    
    if (id<nLifeHistories) {
      outParameters.record("age_d",now());
    }
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
    break;

  case toOrganised:
    organised = true;
    RemoveKind(toScreen); // remove other screens
    scheduleAt(now(), toScreen); // now start organised screening
    break;

  case toBiopsyFollowUpScreen:
  case toScreen:
    if (includePSArecords) {
      psarecord.record("id",id);
      psarecord.record("state",state);
      psarecord.record("ext_grade",ext_grade);
      psarecord.record("organised",organised);
      psarecord.record("dx",dx);
      psarecord.record("age",age);
      psarecord.record("psa",psa);
      psarecord.record("t0",t0);
      psarecord.record("beta0",beta0);
      psarecord.record("beta1",beta1);
      psarecord.record("beta2",beta2);
      psarecord.record("Z",Z);
    }
    costs.add("InvitationCost",now(),cost_parameters["InvitationCost"]);
    // scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["InvitationUtility"]));
    // scheduleAt(now() + utility_duration["InvitationUtilityDuration"], 
    // 	       new cMessageUtilityChange(utility_estimates["InvitationUtility"]));
    if (organised) {
      if (screen == stockholm3_risk_stratified && psa>=1.0) 
	costs.add("FormalPSABiomarkerCost",now(),cost_parameters["FormalPSABiomarkerCost"]);
      else
	costs.add("FormalPSACost",now(),cost_parameters["FormalPSACost"]);
      scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["FormalPSAUtility"]));
      scheduleAt(now() + utility_duration["FormalPSAUtilityDuration"], 
		 new cMessageUtilityChange(utility_estimates["FormalPSAUtility"]));
    } else { // opportunistic
      costs.add("OpportunisticPSACost",now(),cost_parameters["OpportunisticPSACost"]); //cost for opportunistic PSA test
      scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["OpportunisticPSAUtility"]));
      scheduleAt(now() + utility_duration["OpportunisticPSAUtilityDuration"], 
		   new cMessageUtilityChange(utility_estimates["OpportunisticPSAUtility"]));
    }
    if (!everPSA) {
      if (id<nLifeHistories) {
	outParameters.revise("age_psa",now());
      }
      everPSA = true;
    }
    // Other threshold here
    compliance = tableBiopsyCompliance(pair<double,double>(bounds<double>(psa,4.0,7.0),
							   bounds<double>(age,55,75)));
    if (msg->kind == toScreen && psa >= parameter["psaThreshold"] && R::runif(0.0,1.0) < compliance) {
      scheduleAt(now(), toScreenInitiatedBiopsy); // immediate biopsy
    } 
    else if (msg->kind == toBiopsyFollowUpScreen && psa >= parameter["psaThresholdBiopsyFollowUp"] && R::runif(0.0,1.0) < compliance) {
      	scheduleAt(now(), toScreenInitiatedBiopsy); // immediate biopsy, assuming similar biopsy compliance, reasonable? An option to different psa-thresholds would be to use different biopsyCompliance. /AK
    }
    else { // re-screening schedules
      rngScreen->set();
      if (organised) { // NB: currently assumes 100% compliance?
	switch (screen) {
	case stockholm3_goteborg:
	  if (psa<1.0)
	    scheduleAt(now() + 4.0, toScreen);
	  else
	    scheduleAt(now() + 2.0, toScreen);
	  break;
	case stockholm3_risk_stratified:
	  if (psa>=1.0)
	    costs.add("FormalPSABiomarkerCost",now(),cost_parameters["FormalPSABiomarkerCost"]);
	  else
	    costs.add("FormalPSACost",now(),cost_parameters["FormalPSACost"]);
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
      } else  // opportunistic screening - move to a separate function?
	switch(screen) {
	case screenUptake:
	case stockholm3_goteborg:
	case stockholm3_risk_stratified: {
	  TableDDD::key_type key = TableDDD::key_type(bounds<double>(now(),30.0,90.0),psa);
	  double prescreened = 1.0 - rescreen_cure(key); 
	  double shape = rescreen_shape(key);
	  double scale = rescreen_scale(key);
	  double u = R::runif(0.0,1.0);
	  double t = now() + R::rweibull(shape,scale);
	  // Rprintf("(prescreened=%f,u=%f,shape=%f,scale=%f,t=%f)\n",prescreened,u,shape,scale,t);
	  if (u<prescreened) {
	    scheduleAt(t, toScreen);
	  }
	}
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

    if (state == Healthy) {
      previousNegativeBiopsy = true;
      // Here we want 20% to opportunistic and 80% to re-screen in 12 months with threshold of 4. N.B. also the false negative 6 rows below.
      if (now() < 70.0 && R::runif(0.0,1.0)<parameter["screeningCompliance"]) scheduleAt(now() + 1.0, toBiopsyFollowUpScreen);
      // else schedule a routine future screen
    } else { // state != Healthy
      if (state == Metastatic || (state == Localised && R::runif(0.0, 1.0) < parameter["biopsySensitivity"])) {
	scheduleAt(now(), toScreenDiagnosis);
      } else { // false negative biopsy
	if (now() < 70.0 && R::runif(0.0,1.0)<parameter["screeningCompliance"]) scheduleAt(now() + 1.0, toBiopsyFollowUpScreen);
      }
    }
    break;

  case toTreatment: {
    if (state == Metastatic) {
      costs.add("MetastaticCancerCost",now(),cost_parameters["MetastaticCancerCost"]);
      scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["MetastaticCancerUtility"]));
      // scheduleAt(now() + utility_duration["MetastaticCancerUtilityDuration"], 
      // 		 new cMessageUtilityChange(utility_estimates["MetastaticCancerUtility"]));
    }
    else { // Localised
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
    }
    // reset the stream
    rngNh->set();
    // calculate survival
    txhaz = (state == Localised && (tx == RP || tx == RT)) ? 0.62 : 1.0;
    double sxbenefit = 1;
    double u = R::runif(0.0,1.0);
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
    // IHE: Metastatic cancer cf. MISCAN: Palliative treatment
    // Localised -> Metastatic (36/12) -> Cancer death
    if (state == Localised) {
      if (age_cancer_death>now()+utility_duration["MetastaticCancerUtilityDuration"]+utility_duration["PalliativeUtilityDuration"])
	scheduleAt(age_cancer_death - utility_duration["MetastaticCancerUtilityDuration"] - utility_duration["PalliativeUtilityDuration"], 
		   new cMessageUtilityChange(-utility_estimates["MetastaticCancerUtility"]));
      else // cancer death within 36 months of diagnosis/treatment
	scheduleAt(now(),
		   new cMessageUtilityChange(-utility_estimates["MetastaticCancerUtility"]));
    }
    // IHE: Palliative care cf. MISCAN: Terminal care
    // Metastatic -> Palliative (6/12) -> cancer death
    if (age_cancer_death>now()+utility_duration["PalliativeUtilityDuration"])
      scheduleAt(age_cancer_death - utility_duration["PalliativeUtilityDuration"], 
		 new cMessageUtilityChange(-utility_estimates["PalliativeUtility"]+utility_estimates["MetastaticCancerUtility"]));
    else // cancer death within six months of diagnosis/treatment
      scheduleAt(now(),
		 new cMessageUtilityChange(-utility_estimates["PalliativeUtility"]+utility_estimates["MetastaticCancerUtility"]));
      
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
    // costs & utilities??
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
  cost_parameters = as<NumericVector>(otherParameters["cost_parameters"]);
  utility_estimates = as<NumericVector>(otherParameters["utility_estimates"]);
  utility_duration = as<NumericVector>(otherParameters["utility_duration"]);

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
  tableBiopsyCompliance = TableBiopsyCompliance(as<DataFrame>(tables["biopsyComplianceTable"]),
						"psa","age","compliance");
  rescreen_shape = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "shape");
  rescreen_scale = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "scale");
  rescreen_cure  = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "cure");
  
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

  report.discountRate = 0.03;
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
  return List::create(_("costs") = costs.wrap(),                // CostReport
		      _("summary") = report.wrap(),             // EventReport 
		      _("lifeHistories") = wrap(lifeHistories), // vector<LifeHistory::Type>
		      _("parameters") = outParameters.wrap(),   // SimpleReport<double>
		      _("psarecord")=psarecord.wrap()            // SimpleReport<double>
		      );
}

} // anonymous namespace
