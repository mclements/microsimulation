/**
 * @file
 * @author  Mark Clements <mark.clements@ki.se>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * Prostate cancer screening simulation model.
 */


#include "microsimulation.h"

#include <boost/algorithm/cxx11/iota.hpp>

namespace fhcrc_example {

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
		toScreenInitiatedBiopsy,toClinicalDiagnosticBiopsy,toScreenDiagnosis,toOrganised,toTreatment,toCM,toRP,toRT,toADT,toUtilityChange, toUtility, toSTHLM3, toOpportunistic };

  enum screen_t {noScreening, randomScreen50to70, twoYearlyScreen50to70, fourYearlyScreen50to70,
		 screen50, screen60, screen70, screenUptake, stockholm3_goteborg, stockholm3_risk_stratified,
		 goteborg, risk_stratified, mixed_screening,
		 regular_screen, single_screen};

  enum treatment_t {no_treatment, CM, RP, RT};

  enum survival_t {StageShiftBased, LeadTimeBased};

  namespace FullState {
    typedef boost::tuple<short,short,short,bool,double> Type;
    enum Fields {state, ext_grade, dx, psa_ge_3, cohort};
    // string names[5] = {"state","ext_grade","dx","psa_ge_3","cohort"};
  }
  namespace LifeHistory {
    typedef boost::tuple<int,short,short,int,short,double,double,double,double> Type;
    enum Fields {id,state,ext_grade,dx,event,begin,end,year,psa};
  }

  RcppExport SEXP rllogis_(SEXP shape, SEXP scale) {
    RNGScope scope;
    return wrap(R::rllogis(as<double>(shape),as<double>(scale)));
  }
  RcppExport SEXP rllogis_trunc_(SEXP shape, SEXP scale, SEXP left) {
    RNGScope scope;
    return wrap(R::rllogis_trunc(as<double>(shape),as<double>(scale),as<double>(left)));
  }

  EventReport<FullState::Type,short,double> report;
  EventReport<int,short,double> shortReport;
  typedef pair<string,double> CostKey;
  CostReport<CostKey> costs;
  vector<LifeHistory::Type> lifeHistories;
  SimpleReport<double> outParameters;
  SimpleReport<double> psarecord, falsePositives;
  SimpleReport<double> diagnoses;

  bool debug = false;

  typedef Table<boost::tuple<double,double,int>,double> TablePrtx; // Age, DxY, G
  typedef Table<boost::tuple<double,int,int>,double> TableLocoHR; // Age, G, PSA10+
  typedef Table<double,double> TableMetastaticHR; // Age
  typedef Table<boost::tuple<int,double,double,int>,double> TablePradt;
  typedef Table<pair<double,double>,double> TableBiopsyCompliance;
  typedef Table<pair<double,double>,double> TableDDD; // as per TableBiopsyCompliance
  typedef map<int,NumericInterpolate> H_dist_t;
  typedef map<pair<double,int>,NumericInterpolate> H_local_t;
  TableLocoHR hr_locoregional;
  TableMetastaticHR hr_metastatic;
  TablePrtx prtxCM, prtxRP;
  TablePradt pradt;
  TableBiopsyCompliance tableOpportunisticBiopsyCompliance, tableFormalBiopsyCompliance;
  TableDDD rescreen_shape, rescreen_scale, rescreen_cure;
  NumericInterpolate interp_prob_grade7;
  H_dist_t H_dist;
  H_local_t H_local;
  set<double,greater<double> > H_local_age_set;

  Rng * rngNh, * rngOther, * rngScreen, * rngTreatment;
  Rpexp rmu0;

  NumericVector parameter;
  LogicalVector bparameter;

  // read in the parameters
  NumericVector cost_parameters, utility_estimates, utility_duration;
  NumericVector mubeta2, sebeta2; // otherParameters["mubeta2"] rather than as<NumericVector>(otherParameters["mubeta2"])
  int screen, nLifeHistories;
  bool includePSArecords, panel, includeDiagnoses;

  // an alternative approach to specifying costs and compliance (not currently used)
  class cMessageScreen : public cMessage {
  public:
    cMessageScreen(bool formal_compliance = false, bool formal_costs = false) :
      cMessage(toScreen),
      formal_compliance(formal_compliance),
      formal_costs(formal_costs) { }
    bool formal_compliance, formal_costs;
  };

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
    double t0, y0, tm, tc, tmc, aoc;
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
    double psamean(double age);
    double psameasured(double age);
    treatment_t calculate_treatment(double u, double age, double year);
    double calculate_mortality_hr(double age_diag);
    double calculate_survival(double u, double age_diag, double age_c, treatment_t tx);
    void opportunistic_rescreening(double psa);
    void init();
    void add_costs(string item);
    virtual void handleMessage(const cMessage* msg);
    void scheduleUtilityChange(double at, std::string category, bool transient = true,
			       double sign = -1.0);
    bool onset();
  };

  /**
      Calculate the (geometric) mean PSA value at a given age (** NB: this used to be t=age-35.0 **)
  */
  double FhcrcPerson::psamean(double age) {
    double t = age<35.0 ? 0.0 : age - 35.0;
    double yt = t<t0 ? exp(beta0+beta1*t) : exp(beta0+beta1*t+beta2*(t-t0));
    return yt;
  }

  /**
      Calculate the *measured* PSA value at a given age (** NB: this used to be t=age-35 **)
  */
  double FhcrcPerson::psameasured(double age) {
    return FhcrcPerson::psamean(age)*exp(R::rnorm(0.0, sqrt(double(parameter["tau2"]))));
    }

  /**
      Report on costs for a given item
  */
  void FhcrcPerson::add_costs(string item) {
    costs.add(CostKey(item,cohort),now(),cost_parameters[item]);
  }

  /**
     Schedule a transient utility change.
     Default: sign = -1
   **/
  void FhcrcPerson::scheduleUtilityChange(double at, std::string category, bool transient, double sign) {
    scheduleAt(at, new cMessageUtilityChange(sign*utility_estimates[category]));
    if (transient) {
      scheduleAt(at + utility_duration[category],
		 new cMessageUtilityChange(-sign*utility_estimates[category]));
    }
  }

  treatment_t FhcrcPerson::calculate_treatment(double u, double age, double year) {
    TablePrtx::key_type key;
    if (bparameter["stockholmTreatment"])
       key =
	 TablePrtx::key_type(bounds<double>(age,50.0,85.0),
			     bounds<double>(year,2008.0,2012.0),
			     int(ext_grade));
    else // original FHCRC table prtx
      key =
	TablePrtx::key_type(bounds<double>(age,50.0,79.0),
			    bounds<double>(year,1973.0,2004.0),
			    int(grade));
      double pCM = prtxCM(key);
      double pRP = prtxRP(key);
      treatment_t tx = (u<pCM) ? CM :
	           (u<pCM+pRP) ? RP : RT;
      if (debug) Rprintf("id=%i, Age=%3.0f, DxY=%4.0f, stage=%i, grade=%i, tx=%i, u=%8.6f, pCM=%8.6f, pRP=%8.6f\n",id,age,year,state,grade,(int)tx,u,pCM,pRP);
      return tx;
  }

  /** @brief calculate survival taking account of screening and treatment.

      For the Nordic natural history model, we have calibrated
      survival to observed survival from PCBase. These calibrations
      are represented by the hr_locoregional and hr_metastatic tables.

      Note that we do not use now(), as the time points change
      depending on how we represent screening. This explains why we
      pass all of the ages as parameters.

      The ustar is an approach to adjust survival for different
      factors. Rather than solve a revised survival function, we solve
      the baseline survival function for a revised ustar = u^(1/HR).
   **/

  double FhcrcPerson::calculate_mortality_hr(double age_diag) { // also: tm, ext_grade
    double age_m = tm + 35.0;   // age at onset of metastatic cancer
    bool localised = (age_diag < age_m);
    double mort_hr;
    if (localised)
      mort_hr = hr_locoregional(TableLocoHR::key_type(age_diag<50.0 ? 50.0 : age_diag, ext_grade, psamean(age_diag)>10 ? 1 : 0));
    else
      mort_hr = hr_metastatic(age_diag);
    return mort_hr;
  }

  double FhcrcPerson::calculate_survival(double u, double age_diag, double age_c, treatment_t tx) { // also: tc, tm, tmc, grade, now()
    double age_d = -1.0;        // age at death (output)
    double age_m = tm + 35.0;   // age at onset of metastatic cancer
    bool localised = (age_diag < age_m);
    double txhaz = (localised && (tx == RP || tx == RT)) ? 0.62 : 1.0;
    // calibration HR(age_diag,PSA,ext_grade) for loco-regional or HR(age_diag) for metastatic cancer
    double lead_time = age_c - age_diag;
    double txbenefit = exp(log(txhaz)+log(double(parameter["c_txlt_interaction"]))*lead_time);
    double mort_hr = calculate_mortality_hr(age_diag);
    double ustar = pow(u,1/(parameter["c_baseline_specific"]*mort_hr*txbenefit*parameter["sxbenefit"]));
    if (localised)
      age_d = age_c + H_local[H_local_t::key_type(*H_local_age_set.lower_bound(bounds<double>(age_diag,50.0,80.0)),grade)].invert(-log(ustar));
    else
      age_d = age_c + H_dist[grade].invert(-log(ustar));
    if (debug) Rprintf("id=%i, lead_time=%f, ext_grade=%i, psamean=%f, tx=%i, txbenefit=%f, u=%f, ustar=%f, age_diag=%f, age_m=%f, age_c=%f, age_d=%f, mort_hr=%f\n",
		       id, lead_time, (int)ext_grade, psamean(age_diag), (int)tx,txbenefit, u, ustar, age_diag, age_m, age_c, age_d, mort_hr);
    return age_d;
  }

  bool FhcrcPerson::onset() { return now() <= this->t0+35.0; }

  void FhcrcPerson::opportunistic_rescreening(double psa) {
    TableDDD::key_type key = TableDDD::key_type(bounds<double>(now(),30.0,90.0),psa);
    double prescreened = 1.0 - rescreen_cure(key);
    double shape = rescreen_shape(key);
    double scale = rescreen_scale(key);
    double u = R::runif(0.0,1.0);
    double t = now() + R::rweibull(shape,scale);
    if (u<prescreened) {
      scheduleAt(t, toScreen);
    }
  }

  typedef std::pair<double,double> Double;
  Double rbinorm(Double mean, Double sd, double rho) {
    double z1 = R::rnorm(0.0,1.0);
    double z2 = R::rnorm(0.0,1.0);
    z2 = rho*z1+sqrt(1-rho*rho)*z2;
    return Double(z1*sd.first+mean.first, mean.second+sd.second*z2);
  }
  Double rbinormPos(Double mean, Double sd, double rho, Double lbound = Double(0.0,0.0)) {
    Double out;
    do {
      out = rbinorm(mean,sd,rho);
    } while (out.first<lbound.first || out.second<lbound.second);
    return out;
  }
  RcppExport SEXP rbinorm_test() {
    RNGScope rng;
    Double x = rbinorm(Double(1.0,1.0), Double(2.0,2.0), 0.62);
    vector<double> v; v.push_back(x.first); v.push_back(x.second);
    return wrap(v);
  }
  RcppExport SEXP rbinormPos_test() {
    RNGScope rng;
    Double x = rbinormPos(Double(1.0,1.0), Double(2.0,2.0), 0.62);
    vector<double> v; v.push_back(x.first); v.push_back(x.second);
    return wrap(v);
  }


/**
    Initialise a simulation run for an individual
 */
void FhcrcPerson::init() {

  // declarations
  double ym;

  // change state variables
  state = Healthy;
  dx = NotDiagnosed;
  everPSA = previousNegativeBiopsy = organised = adt = false;
  rngNh->set();
  t0 = sqrt(2*R::rexp(1.0)/parameter["g0"]);
  if (!bparameter["revised_natural_history"]){
    grade = (R::runif(0.0, 1.0)>=1+parameter["c_low_grade_slope"]*t0) ? base::Gleason_ge_8 : base::Gleason_le_7;
    beta2 = R::rnormPos(mubeta2[grade],sebeta2[grade]);
  }
  else {
    double u = R::runif(0.0,1.0);
    if (u < exp(parameter["alpha8"] + parameter["beta8"] * t0))
      ext_grade = ext::Gleason_ge_8;
    else if (u > 1 - (parameter["alpha7"] + parameter["beta7"] * t0))
      ext_grade = ext::Gleason_7;
    else ext_grade = ext::Gleason_le_6;
    grade = ext_grade == ext::Gleason_ge_8 ? base::Gleason_ge_8 : base::Gleason_le_7;
    beta2 = R::rnormPos(mubeta2[ext_grade],sebeta2[ext_grade]);
  }
  beta0 = R::rnorm(parameter["mubeta0"],parameter["sebeta0"]);
  beta1 = R::rnormPos(parameter["mubeta1"],parameter["sebeta1"]);

  y0 = psamean(t0+35); // depends on: t0, beta0, beta1, beta2
  tm = (log((beta1+beta2)*R::rexp(1.0)/parameter["gm"] + y0) - beta0 + beta2*t0) / (beta1+beta2);
  ym = psamean(tm+35);
  tc = (log((beta1+beta2)*R::rexp(1.0)/parameter["gc"] + y0) - beta0 + beta2*t0) / (beta1+beta2);
  tmc = (log((beta1+beta2)*R::rexp(1.0)/(parameter["gc"]*parameter["thetac"]) + ym) - beta0 + beta2*t0) / (beta1+beta2);
  aoc = rmu0.rand(R::runif(0.0,1.0));
  if (!bparameter["revised_natural_history"]){
    ext_grade= (grade==base::Gleason_le_7) ?
      (R::runif(0.0,1.0)<=interp_prob_grade7.approx(beta2) ? ext::Gleason_7 : ext::Gleason_le_6) :
      ext::Gleason_ge_8;
  }


  if (debug) {
    Rprintf("id=%i, grade=%i, ext_grade=%i, beta0=%f, beta1=%f, beta2=%f, mubeta0=%f, sebeta0=%f, mubeta1=%f, sebeta1=%f, mubeta2=%f, sebeta2=%f\n", id, grade, ext_grade, beta0, beta1, beta2, double(parameter["mubeta0"]), double(parameter["sebeta0"]), double(parameter["mubeta1"]), double(parameter["sebeta1"]), mubeta2[grade], sebeta2[grade]);
  }

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
    case single_screen:
    case regular_screen:
    case goteborg:
    case risk_stratified:
      scheduleAt(parameter["start_screening"],toScreen);
      break;
    case fourYearlyScreen50to70: // 50,54,58,62,66,70
    case twoYearlyScreen50to70:  // 50,52, ..., 68,70
    case screen50:
      scheduleAt(50.0,toScreen);
      break;
    case screen60:
      scheduleAt(60.0,toScreen);
      break;
    case screen70:
      scheduleAt(70.0,toScreen);
      break;
    case mixed_screening:
    case stockholm3_goteborg:
    case stockholm3_risk_stratified:
    case screenUptake: {
      // Assume:
      // (i)   cohorts aged <35 in 1995 have a llogis(3.8,15) from age 35 (cohort > 1960)
      // (ii)  cohorts aged 50+ in 1995 have a llogis(2,10) distribution from 1995 (cohort < 1945)
      // (iii) intermediate cohorts are a weighted mixture of (i) and (ii)
      double pscreening = cohort>=1932.0 ? 0.9 : 0.9-(1932.0 - cohort)*0.03;
      double shapeA = 3.8;
      double scaleA = 15.0;
      double shapeT = 2.0;
      double scaleT = 10.0;
      double uscreening = R::runif(0.0,1.0);
      double first_screen;
      if (cohort > 1960.0) {
	first_screen = 35.0 + R::rllogis(shapeA,scaleA); // (i) age
      } else if (cohort < 1945.0) {
	first_screen = (1995.0 - cohort) + R::rllogis(shapeT,scaleT); // (ii) period
      } else {
	double age0 = 1995.0 - cohort;
	double u = R::runif(0.0,1.0);
	if ((age0 - 35.0)/15.0 < u) // (iii) mixture
	  first_screen = age0 + R::rllogis_trunc(shapeA,scaleA,age0-35.0);
	else first_screen = age0 + R::rllogis(shapeT,scaleT);
      }
      if (uscreening<pscreening)
	scheduleAt(first_screen, toScreen);
      if (debug)
	Rprintf("(cohort=%f,pscreening=%f,uscreening=%f,first_screen=%f)\n",cohort,pscreening,uscreening,first_screen);
      if (screen == mixed_screening) {
	scheduleAt(parameter["start_screening"], toOrganised);
      }
      // for re-screening patterns: see handleMessage()
    } break;
    default:
      REprintf("Screening not matched: %i\n",screen);
      break;
    }
  }
  if (R::runif(0.0,1.0)<parameter["studyParticipation"] &&
      ((screen == stockholm3_goteborg) || (screen == stockholm3_risk_stratified)) &&
      (2013.0-cohort>=parameter["start_screening"] && 2013.0-cohort<parameter["stop_screening"])) {
    scheduleAt(R::runif(2013.0,2015.0) - cohort, toSTHLM3);
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
    outParameters.record("psa55",psameasured(55.0));
    outParameters.record("psa65",psameasured(65.0));
    outParameters.record("psa75",psameasured(75.0));
    outParameters.record("psa85",psameasured(85.0));
  }
}

/**
    Handle self-messages received
 */
void FhcrcPerson::handleMessage(const cMessage* msg) {

  // by default, use the natural history RNG
  rngNh->set();

  // declarations
  double psa = psameasured(now()); // includes measurement error
  // double test = panel ? biomarker : psa;
  double Z = psamean(now());
  double age = now();
  double year = age + cohort;
  double compliance;
  bool formal_costs = parameter["formal_costs"]==1.0 && (screen != mixed_screening || organised);
  bool formal_compliance = parameter["formal_compliance"]==1.0 && (screen != mixed_screening || organised);

  // record information
  if (parameter["full_report"] == 1.0)
    report.add(FullState::Type(state, ext_grade, dx, psa>=3.0, cohort), msg->kind, previousEventTime, age, utility);
  shortReport.add(1, msg->kind, previousEventTime, age, utility);

  if (id<nLifeHistories) { // only record up to the first n individuals
    lifeHistories.push_back(LifeHistory::Type(id,state,ext_grade,dx,msg->kind,previousEventTime,age,year,psa));
  }

  // handle messages by kind

  switch(msg->kind) {

  case toCancerDeath:
    add_costs("CancerDeath"); // cost for death, should this be zero???
    if (id<nLifeHistories) {
      outParameters.record("age_d",now());
      outParameters.revise("pca_death",1.0);
    }
    Sim::stop_simulation();
    break;

  case toOtherDeath:
    add_costs("Death"); // cost for death, should this be zero???

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
    scheduleAt(tmc+35.0,toClinicalDiagnosis);
    break;

  case toOrganised:
  case toSTHLM3:
    organised = true;
    RemoveKind(toScreen); // remove other screens
    scheduleAt(now(), toScreen); // now start organised screening
    break;

  case toScreen:
  case toBiopsyFollowUpScreen: {
    rngScreen->set();
    if (includePSArecords) {
      psarecord.record("id",id);
      psarecord.record("state",state);
      psarecord.record("ext_grade",ext_grade);
      psarecord.record("organised",organised); // only meaningful for mixed_screening
      psarecord.record("dx",dx);
      psarecord.record("age",age);
      psarecord.record("psa",psa);
      psarecord.record("t0",t0);
      psarecord.record("beta0",beta0);
      psarecord.record("beta1",beta1);
      psarecord.record("beta2",beta2);
      psarecord.record("Z",Z);
    }
    if (!everPSA) {
      if (id<nLifeHistories) {
	outParameters.revise("age_psa",now());
	// outParameters.revise("first_psa",psa);
      }
      everPSA = true;
    }
    if (formal_costs) {
      add_costs("Invitation");
      add_costs(panel && psa>=1.0 ? "FormalPSABiomarker" : "FormalPSA");
      scheduleUtilityChange(now(), "FormalPSA");
    } else { // opportunistic costs
      add_costs(panel && psa>=1.0 ? "OpportunisticPSABiomarker" : "OpportunisticPSA");
      scheduleUtilityChange(now(), "OpportunisticPSA");
    }
    compliance = formal_compliance ?
      tableFormalBiopsyCompliance(pair<double,double>(bounds<double>(psa,4.0,10.0),
						bounds<double>(age,55,75))) :
      tableOpportunisticBiopsyCompliance(pair<double,double>(bounds<double>(psa,4.0,10.0),
						bounds<double>(age,55,75)));
    // bool positive_test =
    //   (!panel && msg->kind == toScreen && psa >= parameter["psaThreshold"]) ? true :
    //   ( panel && msg->kind == toScreen && biomarker >= parameter["BPThreshold"]) ? true :
    //   (!panel && msg->kind == toBiopsyFollowUpScreen && psa >= parameter["psaThresholdBiopsyFollowUp"]) ? true :
    //   ( panel && msg->kind == toBiopsyFollowUpScreen && biomarker >= parameter["BPThresholdBiopsyFollowUp"]) ? true :
    //   false;
    bool positive_test =
      (msg->kind == toScreen && psa >= parameter["psaThreshold"]) ? true :
      (msg->kind == toBiopsyFollowUpScreen && psa >= parameter["psaThresholdBiopsyFollowUp"]) ? true :
      false;
    // Important case: PSA<1 (to check)
    // Reduce false positives wrt Gleason 7+ by 1-rFPF: which BPThreshold?
    if (panel && positive_test) {
      // if (R::runif(0.0,1.0) < 1.0-parameter["rFPF"]) positive_test = false;
      if ((ext_grade == ext::Gleason_le_6 && onset() && psa<parameter["PSA_FP_threshold_GG6"]) // FP GG 6 PSA threshold
	  ||  (!onset() && psa < parameter["PSA_FP_threshold_nCa"])) {// FP no cancer PSA threshold
	positive_test = false; // strong assumption
      }
    }
    if (includePSArecords && !onset() && positive_test) {
      falsePositives.record("id",id);
      falsePositives.record("psa",psa);
      falsePositives.record("age",now());
      falsePositives.record("age0",t0+35.0);
      falsePositives.record("ext_grade",ext_grade);
    }
    // if (panel && !positive_test && t0<now()-35.0 && ext_grade > ext::Gleason_le_6) {
    //   if (R::runif(0.0,1.0) < 1.0-parameter["rTPF"]) positive_test = true;
    // }
    if (positive_test && R::runif(0.0,1.0) < compliance) {
      scheduleAt(now(), toScreenInitiatedBiopsy); // immediate biopsy
    } // assumes similar biopsy compliance, reasonable? An option to different psa-thresholds would be to use different biopsyCompliance. /AK
    else { // re-screening schedules
      if (R::runif(0.0,1.0)<parameter["rescreeningCompliance"]) {
	switch (screen) {
	case mixed_screening:
	case stockholm3_goteborg:
	case goteborg:
	  if (screen != mixed_screening || organised) {
	    if (now() >= parameter["start_screening"] && now()<parameter["stop_screening"]) {
	      if (psa<1.0 && now()+4.0<=parameter["stop_screening"])
		scheduleAt(now() + 4.0, toScreen);
	      else if (psa>=1.0 && now()+2.0<=parameter["stop_screening"])
		scheduleAt(now() + 2.0, toScreen);
	      else if (screen == mixed_screening) {
		organised = false;
		opportunistic_rescreening(psa); // start opportunistic rescreening
	      }
	    }
	  }
	  break;
	case stockholm3_risk_stratified:
	case risk_stratified:
	  if (now() >= parameter["start_screening"]) {
	    if (psa<1.0 && now()+8.0<=parameter["stop_screening"])
	      scheduleAt(now() + 8.0, toScreen);
	    if (psa>=1.0 && now()+4.0<=parameter["stop_screening"])
	      scheduleAt(now() + 4.0, toScreen);
	  }
	  break;
	case regular_screen:
	  if (parameter["start_screening"] <= now() &&
	      now()+parameter["screening_interval"] <= parameter["stop_screening"])
	    scheduleAt(now() + parameter["screening_interval"], toScreen);
	  break;
	case twoYearlyScreen50to70:
	  if (50.0 <= now() && now() < 70.0)
	    scheduleAt(now() + 2.0, toScreen);
	  break;
	case fourYearlyScreen50to70:
	  if (50.0 <= now() && now() < 70.0)
	    scheduleAt(now() + 4.0, toScreen);
	  break;
	case screenUptake:
	case randomScreen50to70:
	case single_screen:
	case screen50:
	case screen60:
	case screen70:
	  break;
	default:
	  REprintf("Screening not matched: %s\n",screen);
	  break;
	}
      } // rescreening compliance
      if (screen == screenUptake || (screen == mixed_screening && !organised))
	opportunistic_rescreening(psa);
    } // rescreening
    rngNh->set();
  } break;

  case toClinicalDiagnosis:
    dx = ClinicalDiagnosis;
    RemoveKind(toMetastatic); // competing events
    RemoveKind(toScreen);
    scheduleAt(now(), toClinicalDiagnosticBiopsy); // for reporting (assumes three biopsies per clinical diagnosis)
    scheduleAt(now(), toClinicalDiagnosticBiopsy);
    scheduleAt(now(), toClinicalDiagnosticBiopsy);
    scheduleAt(now(), toTreatment);
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
    add_costs("Biopsy");
    scheduleUtilityChange(now(), "Biopsy");
    break;

  case toScreenInitiatedBiopsy:
    rngScreen->set();
    add_costs("Biopsy");
    scheduleUtilityChange(now(), "Biopsy");

    if (state == Healthy) {
      previousNegativeBiopsy = true;
      // Here we want 20% to opportunistic and 80% to re-screen in 12 months with threshold of 4. N.B. also the false negative 6 rows below.
      if (now() < parameter["stop_screening"] && R::runif(0.0,1.0)<parameter["screeningCompliance"]) scheduleAt(now() + 1.0, toBiopsyFollowUpScreen);
      // else schedule a routine future screen
    }
    else { // state != Healthy
      if (state == Metastatic || (state == Localised && R::runif(0.0, 1.0) < parameter["biopsySensitivity"])) {
	scheduleAt(now(), toScreenDiagnosis);
      } else { // false negative biopsy
	if (now() < parameter["stop_screening"] && R::runif(0.0,1.0)<parameter["screeningCompliance"]) scheduleAt(now() + 1.0, toBiopsyFollowUpScreen);
      }
    }
    rngNh->set();
    break;

  case toTreatment: {
    rngTreatment->set();
    double u_tx = R::runif(0.0,1.0);
    double u_adt = R::runif(0.0,1.0);
    if (state == Metastatic) {
      RemoveKind(toUtility);
      scheduleAt(now(), new cMessageUtilityChange(-utility_estimates["MetastaticCancer"]));
    }
    else { // Loco-regional
      tx = calculate_treatment(u_tx,now(),year);
      if (tx == CM) scheduleAt(now(), toCM);
      if (tx == RP) scheduleAt(now(), toRP);
      if (tx == RT) scheduleAt(now(), toRT);
      // check for ADT
      double pADT =
	pradt(TablePradt::key_type(tx,
				   bounds<double>(now(),50,79),
				   bounds<double>(year,1973,2004),
				   grade));
      if (u_adt < pADT)  {
	adt = true;
	scheduleAt(now(), toADT);
      }
      if (debug) Rprintf("id=%i, adt=%d, u=%8.6f, pADT=%8.6f\n",id,adt,u_adt,pADT);
    }
    // reset the random number stream
    rngNh->set();
    // check for cure
    bool cured = false;
    double age_c = (state == Localised) ? tc + 35.0 : tmc + 35.0;
    double lead_time = age_c - now();
    // calculate the age at cancer death by c_benefit_type
    double age_cancer_death=R_PosInf;
    if (parameter["c_benefit_type"]==LeadTimeBased) { // [new paper ref]
      double pcure = pow(1 - exp(-lead_time*parameter["c_benefit_value1"]),
      			 calculate_mortality_hr(age_c));
      if (debug) Rprintf("hr for lead time=%f\n", calculate_mortality_hr(age_c));
      cured = (R::runif(0.0,1.0) < pcure);
      if (cured) RemoveKind(toMetastatic);
      else {
      double u_surv = R::runif(0.0,1.0);
      age_cancer_death = calculate_survival(u_surv,age_c,age_c,calculate_treatment(u_tx,age_c,year+lead_time));
      }
    }
    else if (parameter["c_benefit_type"]==StageShiftBased) { // [annals paper ref]
      // calculate survival
      double u_surv = R::runif(0.0,1.0);
      double age_cd = calculate_survival(u_surv,age_c,age_c,calculate_treatment(u_tx,age_c,year+lead_time));
      double age_sd = calculate_survival(u_surv,now(),age_c,tx);
      double weight = exp(-parameter["c_benefit_value0"]*lead_time);
      age_cancer_death = weight*age_cd + (1.0-weight)*age_sd;
    }
    else REprintf("c_benefit_type not matched.");
    if (!cured) {
      scheduleAt(age_cancer_death, toCancerDeath);
      // Disutilities prior to a cancer death
      double age_palliative = age_cancer_death - utility_duration["PalliativeTherapy"] - utility_duration["TerminalIllness"];
      double age_terminal = age_cancer_death - utility_duration["TerminalIllness"];
      // Reset utilities for those in with a Metatatic diagnosis
      if (state == Metastatic) {
	if (age_palliative > now())
	  scheduleUtilityChange(age_palliative, "MetastaticCancer",false);
	else
	  scheduleUtilityChange(now(), "MetastaticCancer", false);
      }
      if (age_palliative>now()) { // cancer death more than 36 months after diagnosis
	scheduleUtilityChange(age_palliative, "PalliativeTherapy");
	scheduleUtilityChange(age_terminal, "TerminalIllness");
      }
      else if (age_terminal>now()) { // cancer death between 36 and 6 months of diagnosis
	scheduleUtilityChange(now(), "PalliativeTherapy",false, -1.0);
	scheduleUtilityChange(age_terminal, "PalliativeTherapy", false, 1.0); // reset
	scheduleUtilityChange(age_terminal,"TerminalIllness");
      }
      else // cancer death within 6 months of diagnosis/treatment
	scheduleUtilityChange(now(), "TerminalIllness");
    }
    if (includeDiagnoses) {
      diagnoses.record("id",id);
      diagnoses.record("age",age);
      diagnoses.record("year",year);
      diagnoses.record("psa",psa);
      diagnoses.record("ext_grade",ext_grade);
      diagnoses.record("state",state);
      diagnoses.record("organised",organised); // only meaningful for mixed_screening, keep this?
      diagnoses.record("dx",dx);
      diagnoses.record("tx",tx);
      diagnoses.record("cancer_death",(aoc>age_cancer_death) ? 1.0 : 0.0);
      diagnoses.record("age_at_death", (aoc>age_cancer_death) ? age_cancer_death : aoc);
    }
  } break;

  case toRP:
    add_costs("Prostatectomy");
    // Scheduling utilities for the first 2 months after procedure
    scheduleUtilityChange(now(), "ProstatectomyPart1");
    // Scheduling utilities for the first 3-12 months after procedure
    scheduleUtilityChange(now() + utility_duration["ProstatectomyPart1"],
			  "ProstatectomyPart2");
    break;

  case toRT:
    add_costs("RadiationTherapy");
    // Scheduling utilities for the first 2 months after procedure
    scheduleUtilityChange(now(), "RadiationTherapyPart1");
    // Scheduling utilities for the first 3-12 months after procedure
    scheduleUtilityChange(now() + utility_duration["RadiationTherapyPart1"],
			  "RadiationTherapyPart2");
    break;

  case toCM:
    add_costs("ActiveSurveillance"); // expand here
    scheduleUtilityChange(now(), "ActiveSurveillance");
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
  bparameter = parms["bparameter"]; // scalar bools
  List otherParameters = parms["otherParameters"];
  debug = as<bool>(parms["debug"]);
  if (!bparameter["revised_natural_history"]) {
    mubeta2 = as<NumericVector>(otherParameters["mubeta2"]);
    sebeta2 = as<NumericVector>(otherParameters["sebeta2"]);
  } else {
    mubeta2 = as<NumericVector>(otherParameters["rev_mubeta2"]);
    sebeta2 = as<NumericVector>(otherParameters["rev_sebeta2"]);
  }
  NumericVector mu0 = as<NumericVector>(otherParameters["mu0"]);
  cost_parameters = as<NumericVector>(otherParameters["cost_parameters"]);
  utility_estimates = as<NumericVector>(otherParameters["utility_estimates"]);
  utility_duration = as<NumericVector>(otherParameters["utility_duration"]);

  int n = as<int>(parms["n"]);
  includePSArecords = as<bool>(parms["includePSArecords"]);
  includeDiagnoses = as<bool>(parms["includeDiagnoses"]);
  int firstId = as<int>(parms["firstId"]);
  interp_prob_grade7 =
    NumericInterpolate(as<DataFrame>(tables["prob_grade7"]));
  prtxCM = TablePrtx(as<DataFrame>(tables["prtx"]),
			       "Age","DxY","G","CM"); // NB: Grade is now {0,1[,2]} coded cf {1,2[,3]}
  prtxRP = TablePrtx(as<DataFrame>(tables["prtx"]),
			       "Age","DxY","G","RP");
  pradt = TablePradt(as<DataFrame>(tables["pradt"]),"Tx","Age","DxY","Grade","ADT");
  hr_locoregional = TableLocoHR(as<DataFrame>(otherParameters["hr_locoregional"]),"age","ext_grade","psa10","hr");
  hr_metastatic = TableMetastaticHR(as<DataFrame>(otherParameters["hr_metastatic"]),"age","hr");
  tableOpportunisticBiopsyCompliance = TableBiopsyCompliance(as<DataFrame>(tables["biopsyOpportunisticComplianceTable"]),
						"psa","age","compliance");
  tableFormalBiopsyCompliance = TableBiopsyCompliance(as<DataFrame>(tables["biopsyFormalComplianceTable"]),
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
    // Rprintf("Biopsy compliance: %f\n",tableBiopsyCompliance(pair<double,double>(bounds<double>(1.0,4.0,10.0), bounds<double>(100.0,55,75))));
    Rprintf("Interp for grade 6/7 (expecting approx 0.3): %f\n",interp_prob_grade7.approx(0.143));
    Rprintf("prtxCM(80,2008,1) [expecting 0.970711]: %f\n",prtxCM(TablePrtx::key_type(80.0,2008.0,1)));
    {
      double age_diag=51.0;
      ext::grade_t ext_grade = ext::Gleason_ge_8;
      // FhrcPerson person = FhcrcPerson(0,1960);
      // Rprintf("hr_localregional(50,0,)=%g\n",hr_locoregional(TableLocoHR::key_type(age_diag<50.0 ? 50.0 : age_diag, ext_grade, person.psamean(age_diag)>10 ? 1 : 0)));
      Rprintf("hr_localregional(50,8+,0)=%g\n",hr_locoregional(TableLocoHR::key_type(age_diag<50.0 ? 50.0 : age_diag, ext_grade, 0)));
      Rprintf("hr_localregional(50,8+,1)=%g\n",hr_locoregional(TableLocoHR::key_type(age_diag<50.0 ? 50.0 : age_diag, ext_grade, 1)));
      Rprintf("hr_localregional(50,7,0)=%g\n",hr_locoregional(TableLocoHR::key_type(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_7, 0)));
      Rprintf("hr_localregional(50,7,1)=%g\n",hr_locoregional(TableLocoHR::key_type(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_7, 1)));
      Rprintf("hr_localregional(50,<=6,0)=%g\n",hr_locoregional(TableLocoHR::key_type(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_le_6, 0)));
      Rprintf("hr_localregional(50,<=6,1)=%g\n",hr_locoregional(TableLocoHR::key_type(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_le_6, 1)));
    }
  }

  nLifeHistories = as<int>(otherParameters["nLifeHistories"]);
  screen = as<int>(otherParameters["screen"]);
  panel = as<bool>(parms["panel"]);
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
  shortReport.clear();
  costs.clear();
  outParameters.clear();
  lifeHistories.clear();
  psarecord.clear();
  falsePositives.clear();
  diagnoses.clear();

  report.discountRate = parameter["discountRate.effectiveness"];
  report.setPartition(ages);
  shortReport.discountRate = parameter["discountRate.effectiveness"];
  shortReport.setPartition(ages);
  costs.discountRate = parameter["discountRate.costs"];
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
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
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
		      _("shortSummary") = shortReport.wrap(),   // EventReport
		      _("lifeHistories") = wrap(lifeHistories), // vector<LifeHistory::Type>
		      _("parameters") = outParameters.wrap(),   // SimpleReport<double>
		      _("psarecord")=psarecord.wrap(),          // SimpleReport<double>
		      _("falsePositives")=falsePositives.wrap(),// SimpleReport<double>
		      _("diagnoses")=diagnoses.wrap()           // SimpleReport<double>
		      );
}

} // anonymous namespace
