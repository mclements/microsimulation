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
    enum grade_t {Gleason_le_7,Gleason_ge_8,Healthy};
  }
  namespace ext {
    enum grade_t {Gleason_le_6,Gleason_7,Gleason_ge_8,Healthy};
    enum state_t {Healthy_state,T1_T2,T3plus,Metastatic};
  }

  enum state_t {Healthy,Localised,Metastatic}; // stage?

  enum diagnosis_t {NotDiagnosed,ClinicalDiagnosis,ScreenDiagnosis};

  enum event_t {toLocalised,toMetastatic,toClinicalDiagnosis,toCancerDeath,toOtherDeath,toScreen,toBiopsyFollowUpScreen,
		toScreenInitiatedBiopsy,toClinicalDiagnosticBiopsy,toScreenDiagnosis,toOrganised,toTreatment,toCM,toRP,toRT,toADT,toUtilityChange, toBaselineUtility, toSTHLM3, toOpportunistic, toT3plus };

  enum screen_t {noScreening, randomScreen50to70, twoYearlyScreen50to70, fourYearlyScreen50to70,
		 screen50, screen60, screen70, screenUptake, stockholm3_goteborg, stockholm3_risk_stratified,
		 goteborg, risk_stratified, mixed_screening,
		 regular_screen, single_screen};

  enum treatment_t {no_treatment, CM, RP, RT};

  enum survival_t {StageShiftBased, LeadTimeBased};

  enum biomarker_model_t {random_correction, psa_informed_correction};

  enum cost_t {Direct,Indirect};

  namespace FullState {
    typedef boost::tuple<short,short,short,bool,double> Type;
    enum Fields {ext_state, ext_grade, dx, psa_ge_3, cohort};
    // string names[5] = {"ext_state","ext_grade","dx","psa_ge_3","cohort"};
  }
  namespace LifeHistory {
    typedef boost::tuple<int, short, short, int, short, double, double, double, double, double> Type;
    enum Fields {id, ext_state, ext_grade, dx, event, begin, end, year, psa, utility};
  }

  RcppExport SEXP rllogis_(SEXP shape, SEXP scale) {
    RNGScope scope;
    return wrap(R::rllogis(as<double>(shape),as<double>(scale)));
  }
  RcppExport SEXP rllogis_trunc_(SEXP shape, SEXP scale, SEXP left) {
    RNGScope scope;
    return wrap(R::rllogis_trunc(as<double>(shape),as<double>(scale),as<double>(left)));
  }

  typedef std::pair<int,string> CostKey; // (cost_type,cost_name)
  // typedef boost::tuple<int,string,double> CostKey; // (cost_type,cost_name,cohort)
  class SimOutput {
  public:
    EventReport<FullState::Type,short,double> report;
    EventReport<int,short,double> shortReport;
    CostReport<CostKey> costs;
    vector<LifeHistory::Type> lifeHistories;
    SimpleReport<double> outParameters;
    SimpleReport<double> psarecord, falsePositives;
    SimpleReport<double> diagnoses;
    Means tmc_minus_t0;
  };
  // SimOutput * out; // in callFhcrc
  // &out in the Person object
  // out->report etc


  typedef std::pair<double,double> Double;
  typedef Table<double,double,int,double> TablePrtx; // Age, DxY, G
  typedef Table<double,int,int,double> TableLocoHR; // Age, G, PSA10+
  typedef Table<double,double> TableMetastaticHR; // Age
  typedef Table<int,double,double,int,double> TablePradt;
  typedef Table<double,double,double> TableBiopsyCompliance;
  typedef Table<double,double,double> TableDDD; // as per TableBiopsyCompliance
  typedef map<int,NumericInterpolate> H_dist_t;
  typedef map<pair<double,int>,NumericInterpolate> H_local_t;
  class SimInput {
  public:
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
    NumericVector cost_parameters, utility_estimates, utility_duration, lost_production_proportions;
    NumericVector mubeta2, sebeta2; // otherParameters["mubeta2"] rather than as<NumericVector>(otherParameters["mubeta2"])
    int screen, nLifeHistories;
    bool includePSArecords, panel, includeDiagnoses, debug;
    Table<double,double> production;

    ~SimInput() {
      if (rngNh != NULL) delete rngNh;
      if (rngOther != NULL) delete rngOther;
      if (rngScreen != NULL) delete rngScreen;
      if (rngTreatment != NULL) delete rngTreatment;
    }
  };
  // SimInput in; // callFhcrc
  // SimInput* in; // Person
  // in->hr_locoregional(...);

  class cMessageUtilityChange : public cMessage {
  public:
    cMessageUtilityChange(double change) : cMessage(toUtilityChange), change(change) { }
    double change;
  };

  class cMessageBaselineUtility : public cMessage {
  public:
    cMessageBaselineUtility(double utility) : cMessage(toBaselineUtility), utility(utility) { }
    double utility;
  };

  template<class T>
  T bounds(T x, T a, T b) {
    return (x<a)?a:((x>b)?b:x);
  }

  class FhcrcPerson : public cProcess
  {
  public:
    SimInput* in;
    SimOutput* out;
    double beta0, beta1, beta2;
    double t0, y0, t3p, tm, tc, tmc, aoc;
    state_t state;
    ext::state_t ext_state;
    diagnosis_t dx;
    base::grade_t future_grade, grade;
    ext::grade_t future_ext_grade, ext_grade;
    treatment_t tx;
    bool adt;
    double txhaz;
    int id;
    double cohort, baseline_utility, delta_utility;
    bool everPSA, previousNegativeBiopsy, organised;
    FhcrcPerson(SimInput* in, SimOutput* out, const int id = 0, const double cohort = 1950) :
      in(in), out(out), id(id), cohort(cohort), baseline_utility(1.0), delta_utility(0.0) { };
    double utility() { return baseline_utility + delta_utility; }
    double psamean(double age);
    double psameasured(double age);
    treatment_t calculate_treatment(double u, double age, double year);
    double calculate_mortality_hr(double age_diag);
    double calculate_survival(double u, double age_diag, double age_c, treatment_t tx);
    double calculate_transition_time(double u, double t_enter, double gamma);
    void opportunistic_rescreening(double psa);
    void opportunistic_uptake();
    void init();
    void add_costs(string item, cost_t cost_type = Direct);
    void lost_productivity(string item);
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
    return FhcrcPerson::psamean(age)*exp(R::rnorm(0.0, sqrt(double(in->parameter["tau2"]))));
    }

  /**
      Report on costs for a given item
  */
  void FhcrcPerson::add_costs(string item, cost_t cost_type) {
    out->costs.add(CostKey((int) cost_type,item),now(),in->cost_parameters[item]);
  }

  /**
      Report on lost productivity
  */
  void FhcrcPerson::lost_productivity(string item) {
    double loss = in->lost_production_proportions[item] * in->production(now());
    out->costs.add(CostKey((int) Indirect,item),now(),loss);
  }

  /**
     Schedule a transient utility change.
     Default: sign = -1
   **/
  void FhcrcPerson::scheduleUtilityChange(double at, std::string category, bool transient, double sign) {
    scheduleAt(at, new cMessageUtilityChange(sign * in->utility_estimates[category]));
    if (transient) {
      scheduleAt(at + in->utility_duration[category],
		 new cMessageUtilityChange(-sign * in->utility_estimates[category]));
    }
  }

  treatment_t FhcrcPerson::calculate_treatment(double u, double age, double year) {
    double pCM, pRP;
    if (in->bparameter["stockholmTreatment"]) {
       pCM = in->prtxCM(bounds<double>(age,50.0,85.0),
		    bounds<double>(year,2008.0,2012.0),
		    int(ext_grade));
       pRP = in->prtxRP(bounds<double>(age,50.0,85.0),
		    bounds<double>(year,2008.0,2012.0),
		    int(ext_grade));
    } else { // original FHCRC table prtx
       pCM = in->prtxCM(bounds<double>(age,50.0,79.0),
		    bounds<double>(year,1973.0,2004.0),
		    int(grade));
       pRP = in->prtxRP(bounds<double>(age,50.0,79.0),
		    bounds<double>(year,1973.0,2004.0),
		    int(grade));
    }
    treatment_t tx = (u<pCM) ? CM :
      (u<pCM+pRP) ? RP : RT;
    if (in->debug) Rprintf("id=%i, Age=%3.0f, DxY=%4.0f, stage=%i, grade=%i, tx=%i, u=%8.6f, pCM=%8.6f, pRP=%8.6f\n",id,age,year,state,grade,(int)tx,u,pCM,pRP);
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
    if (localised) {
      mort_hr = in->hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext_grade, psamean(age_diag)>10 ? 1 : 0);
      if (ext_state == ext::T3plus) mort_hr*=double(in->parameter["RR_T3plus"]);
    }
    else
      mort_hr = in->hr_metastatic(age_diag);
    return mort_hr;
  }

  double FhcrcPerson::calculate_survival(double u, double age_diag, double age_c, treatment_t tx) { // also: tc, tm, tmc, grade, now()
    double age_d = -1.0;        // age at death (output)
    double age_m = tm + 35.0;   // age at onset of metastatic cancer
    bool localised = (age_diag < age_m);
    double txhaz = (localised && (tx == RP || tx == RT)) ? 0.62 : 1.0; // hazard ratio from SPCG-4
    // calibration HR(age_diag,PSA,ext_grade) for loco-regional or HR(age_diag) for metastatic cancer
    double lead_time = age_c - age_diag;
    double txbenefit = exp(log(txhaz)+log(double(in->parameter["c_txlt_interaction"]))*lead_time); // treatment lead-time interaction
    double mort_hr = calculate_mortality_hr(age_diag);
    double ustar = pow(u,1/(in->parameter["c_baseline_specific"]*mort_hr*txbenefit*in->parameter["sxbenefit"]));
    if (localised)
      age_d = age_c + in->H_local[H_local_t::key_type(* in->H_local_age_set.lower_bound(bounds<double>(age_diag,50.0,80.0)),grade)].invert(-log(ustar));
    else
      age_d = age_c + in->H_dist[grade].invert(-log(ustar));
    if (in->debug) Rprintf("id=%i, lead_time=%f, ext_grade=%i, psamean=%f, tx=%i, txbenefit=%f, u=%f, ustar=%f, age_diag=%f, age_m=%f, age_c=%f, age_d=%f, mort_hr=%f\n",
		       id, lead_time, (int)ext_grade, psamean(age_diag), (int)tx, txbenefit, u, ustar, age_diag, age_m, age_c, age_d, mort_hr);
    return age_d;
  }

  bool FhcrcPerson::onset() { return now() <= this->t0+35.0; }

  /** @brief Calculate transition times for h(t) = y(t)*gamma = exp(beta0+beta1*t+beta2*(t-t0))*gamma
      This is equivalent to solving H(t) = gamma/(beta1+beta2)*(y(t)-y(s)) = -log(U) for entry time s.
  **/
  double FhcrcPerson::calculate_transition_time(double u, double t_enter, double gamma) {
    double y_enter = psamean(35.0 + t_enter);
    return (log(-log(u)*(beta1+beta2)/gamma + y_enter) - beta0 + beta2*t0) / (beta1+beta2);
  }
  
  void FhcrcPerson::opportunistic_rescreening(double psa) {
    double prescreened = 1.0 - in->rescreen_cure(bounds<double>(now(),30.0,90.0),psa);
    double shape = in->rescreen_shape(bounds<double>(now(),30.0,90.0),psa);
    double scale = in->rescreen_scale(bounds<double>(now(),30.0,90.0),psa);
    double u = R::runif(0.0,1.0);
    double t = now() + R::rweibull(shape,scale);
    if (u<prescreened) {
      scheduleAt(t, toScreen);
    }
  }

  void FhcrcPerson::opportunistic_uptake() {
    // Assume:
    // (i)   cohorts aged <35 in 1995 have a llogis(3.8,15) from age 35 (cohort > 1960)
    // (ii)  cohorts aged 50+ in 1995 have a llogis(2,10) distribution from 1995 (cohort < 1945)
    // (iii) intermediate cohorts are a weighted mixture of (i) and (ii)
    double pscreening = cohort>=1932.0 ? 0.9 : 0.9-(1932.0 - cohort)*0.03;
    double shapeA = 3.8;
    double scaleA = 15.0;
    double shapeT = 2.16;
    double scaleT = 11.7;
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
    if (in->debug)
      Rprintf("(cohort=%f,pscreening=%f,uscreening=%f,first_screen=%f)\n",cohort,pscreening,uscreening,first_screen);
  }

  Double rbinorm(Double mean, Double sd, double rho) {
    double z1 = R::rnorm(0.0,1.0);
    double z2 = R::rnorm(0.0,1.0);
    z2 = rho*z1+sqrt(1-rho*rho)*z2;
    return Double(z1*sd.first+mean.first, mean.second+sd.second*z2);
  }
  Double rbinormPos(Double mean, Double sd, double rho, Double lbound = Double(0.0,0.0)) {
    Double val;
    do {
      val = rbinorm(mean,sd,rho);
    } while (val.first<lbound.first || val.second<lbound.second);
    return val;
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
  ext_state = ext::Healthy_state;
  grade = base::Healthy;
  ext_grade = ext::Healthy;
  dx = NotDiagnosed;
  everPSA = previousNegativeBiopsy = organised = adt = false;
  in->rngNh->set();
  t0 = sqrt(2*R::rexp(1.0)/in->parameter["g0"]);
  if (!in->bparameter["revised_natural_history"]){
    future_grade = (R::runif(0.0, 1.0)>=1+in->parameter["c_low_grade_slope"]*t0) ? base::Gleason_ge_8 : base::Gleason_le_7;
    beta2 = R::rnormPos(in->mubeta2[future_grade],in->sebeta2[future_grade]);
  }
  else {
    double u = R::runif(0.0,1.0);
    if (u < exp(in->parameter["alpha8"] + in->parameter["beta8"] * t0))
      future_ext_grade = ext::Gleason_ge_8;
    else if (u > 1 - (in->parameter["alpha7"] + in->parameter["beta7"] * t0))
      future_ext_grade = ext::Gleason_7;
    else future_ext_grade = ext::Gleason_le_6;
    future_grade = future_ext_grade == ext::Gleason_ge_8 ? base::Gleason_ge_8 : base::Gleason_le_7;
    beta2 = R::rnormPos(in->mubeta2[future_ext_grade],in->sebeta2[future_ext_grade]);
  }
  beta0 = R::rnorm(in->parameter["mubeta0"],in->parameter["sebeta0"]);
  beta1 = R::rnormPos(in->parameter["mubeta1"],in->parameter["sebeta1"]);

  y0 = psamean(t0+35); // depends on: t0, beta0, beta1, beta2
  t3p = calculate_transition_time(R::runif(0.0,1.0), t0, in->parameter["g3p"]);
  tm = calculate_transition_time(R::runif(0.0,1.0), t3p, in->parameter["gm"]);
  ym = psamean(tm+35);
  if (future_grade==base::Gleason_le_7) { // Annals
    tc = calculate_transition_time(R::runif(0.0,1.0), t0, in->parameter["gc"]);
    tmc = calculate_transition_time(R::runif(0.0,1.0), tm, in->parameter["gc"]*in->parameter["thetac"]);
  } else {
    tc = calculate_transition_time(R::runif(0.0,1.0), t0, in->parameter["gc"]*in->parameter["grade.clinical.rate.high"]);
    tmc = calculate_transition_time(R::runif(0.0,1.0), tm, in->parameter["gc"]*in->parameter["thetac"]*in->parameter["grade.clinical.rate.high"]);
  }
  aoc = in->rmu0.rand(R::runif(0.0,1.0));
  if (!in->bparameter["revised_natural_history"]){
    future_ext_grade= (future_grade==base::Gleason_le_7) ?
      (R::runif(0.0,1.0) <= in->interp_prob_grade7.approx(beta2) ? ext::Gleason_7 : ext::Gleason_le_6) :
      ext::Gleason_ge_8;
  }

  if (in->debug) {
    Rprintf("id=%i, future_grade=%i, future_ext_grade=%i, beta0=%f, beta1=%f, beta2=%f, mubeta0=%f, sebeta0=%f, mubeta1=%f, sebeta1=%f, mubeta2=%f, sebeta2=%f\n", id, future_grade, future_ext_grade, beta0, beta1, beta2, double(in->parameter["mubeta0"]), double(in->parameter["sebeta0"]), double(in->parameter["mubeta1"]), double(in->parameter["sebeta1"]), in->mubeta2[future_grade], in->sebeta2[future_grade]);
  }

  tx = no_treatment;
  txhaz = -1.0;
  // schedule natural history events
  scheduleAt(t0+35.0,toLocalised);
  scheduleAt(aoc,toOtherDeath);

  // schedule screening events that depend on screeningCompliance
  in->rngScreen->set();
  if (R::runif(0.0,1.0)<in->parameter["screeningCompliance"]) {
    switch(in->screen) {
    case noScreening:
      break; // no screening
    case randomScreen50to70:
      scheduleAt(R::runif(50.0,70.0),toScreen);
      break;
    case single_screen:
    case regular_screen:
    case goteborg:
    case risk_stratified:
      scheduleAt(in->parameter["start_screening"],toScreen);
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
    case screenUptake:
      // see below (models compliance)
      break;
    default:
      REprintf("Screening not matched: %i\n",in->screen);
      break;
    }
  }
  //
  // schedule screening events that already incorporate screening compliance
  switch(in->screen) {
  case mixed_screening:
    opportunistic_uptake();
    scheduleAt(in->parameter["start_screening"], toOrganised);
    break;
  case stockholm3_goteborg:
  case stockholm3_risk_stratified:
    opportunistic_uptake();
    if (R::runif(0.0,1.0) < in->parameter["studyParticipation"] &&
	(2013.0-cohort >= in->parameter["start_screening"] && 2013.0-cohort < in->parameter["stop_screening"]))
      scheduleAt(R::runif(2013.0,2015.0) - cohort, toSTHLM3);
    break;
  case screenUptake:
    opportunistic_uptake();
    break;
  default:
    break;
  }

  in->rngNh->set();

  // utilities
  // | LowerAge | UpperAge | Males | Females |
  // |----------+----------+-------+---------|
  // |        0 |       17 |  1.00 |    1.00 |
  // |       18 |       24 |  0.89 |    0.83 |
  // |       25 |       29 |  0.89 |    0.84 |
  // |       30 |       34 |  0.88 |    0.85 |
  // |       35 |       39 |  0.87 |    0.83 |
  // |       40 |       44 |  0.84 |    0.81 |
  // |       45 |       49 |  0.84 |    0.80 |
  // |       50 |       54 |  0.83 |    0.78 |
  // |       55 |       59 |  0.83 |    0.77 |
  // |       60 |       64 |  0.82 |    0.77 |
  // |       65 |       69 |  0.83 |    0.79 |
  // |       70 |       74 |  0.81 |    0.75 |
  // |       75 |       79 |  0.79 |    0.73 |
  // |       80 |       84 |  0.74 |    0.69 |
  // (i) set initial baseline utility
  baseline_utility = 1.00;
  // (ii) schedule changes in the baseline utility by age
  // Burström and Rehnberg (2006)
  // (require 'cl)
  // (let ((utilities
  // 	 (list 1 0.89 0.89 0.88 0.87 0.84 0.84 0.83 0.83 0.82 0.83 0.81 0.79 0.74))
  // 	(ages
  // 	 (append (list 0 18) (loop for i from 25 to 80 by 5 collect i))))
  //  (loop for utility in utilities for age in ages
  //   do (message (format "scheduleAt(%g, new cMessageBaselineUtility(%g));" age utility))))
  scheduleAt(0, new cMessageBaselineUtility(1));
  scheduleAt(18, new cMessageBaselineUtility(0.89));
  scheduleAt(25, new cMessageBaselineUtility(0.89));
  scheduleAt(30, new cMessageBaselineUtility(0.88));
  scheduleAt(35, new cMessageBaselineUtility(0.87));
  scheduleAt(40, new cMessageBaselineUtility(0.84));
  scheduleAt(45, new cMessageBaselineUtility(0.84));
  scheduleAt(50, new cMessageBaselineUtility(0.83));
  scheduleAt(55, new cMessageBaselineUtility(0.83));
  scheduleAt(60, new cMessageBaselineUtility(0.82));
  scheduleAt(65, new cMessageBaselineUtility(0.83));
  scheduleAt(70, new cMessageBaselineUtility(0.81));
  scheduleAt(75, new cMessageBaselineUtility(0.79));
  scheduleAt(80, new cMessageBaselineUtility(0.74));

  // record some parameters using SimpleReport - too many for a tuple
  if (id < in->nLifeHistories) {
    out->outParameters.record("id",double(id));
    out->outParameters.record("beta0",beta0);
    out->outParameters.record("beta1",beta1);
    out->outParameters.record("beta2",beta2);
    out->outParameters.record("t0",t0);
    out->outParameters.record("tm",tm);
    out->outParameters.record("tc",tc);
    out->outParameters.record("tmc",tmc);
    out->outParameters.record("y0",y0);
    out->outParameters.record("ym",ym);
    out->outParameters.record("aoc",aoc);
    out->outParameters.record("cohort",cohort);
    out->outParameters.record("future_ext_grade",future_ext_grade);
    out->outParameters.record("ext_grade",ext_grade);
    out->outParameters.record("age_psa",-1.0);
    out->outParameters.record("pca_death",0.0);
    out->outParameters.record("psa55",psameasured(55.0));
    out->outParameters.record("psa65",psameasured(65.0));
    out->outParameters.record("psa75",psameasured(75.0));
    out->outParameters.record("psa85",psameasured(85.0));
  }

  if (in->debug) Rprint_actions();

}

/**
    Handle self-messages received
 */
void FhcrcPerson::handleMessage(const cMessage* msg) {

  // by default, use the natural history RNG
  in->rngNh->set();

  // declarations
  double psa = psameasured(now()); // includes measurement error
  // double test = panel ? biomarker : psa;
  double Z = psamean(now());
  double age = now();
  double year = age + cohort;
  double compliance;
  bool formal_costs = in->parameter["formal_costs"]==1.0 && (in->screen != mixed_screening || organised);
  bool formal_compliance = in->parameter["formal_compliance"]==1.0 && (in->screen != mixed_screening || organised);

  // record information
  if (in->parameter["full_report"] == 1.0)
    out->report.add(FullState::Type(ext_state, ext_grade, dx, psa>=3.0, cohort), msg->kind, previousEventTime, age, utility());
  out->shortReport.add(1, msg->kind, previousEventTime, age, utility());

  if (id < in->nLifeHistories) { // only record up to the first n individuals
    out->lifeHistories.push_back(LifeHistory::Type(id, ext_state, ext_grade, dx, msg->kind, previousEventTime, age, year, psa, utility()));
  }

  // handle messages by kind

  switch(msg->kind) {

  case toCancerDeath:
    add_costs("Cancer death"); // cost for death, should this be zero???
    if (id < in->nLifeHistories) {
      out->outParameters.record("age_d",now());
      out->outParameters.revise("pca_death",1.0);
    }
    Sim::stop_simulation();
    break;

  case toOtherDeath:
    // add_costs("Death"); // cost for death, should this be zero???

    if (id < in->nLifeHistories) {
      out->outParameters.record("age_d",now());
    }
    Sim::stop_simulation();
    break;

  case toLocalised:
    state = Localised; ext_state = ext::T1_T2;
    ext_grade = future_ext_grade;
    grade = future_grade;
    scheduleAt(tc+35.0,toClinicalDiagnosis);
    scheduleAt(tm+35.0,toMetastatic);
    scheduleAt(t3p+35.0,toT3plus);
    break;

  case toT3plus:
    ext_state = ext::T3plus;
    break;

  case toMetastatic:
    state = Metastatic; ext_state = ext::Metastatic;
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
    in->rngScreen->set();
    if (in->includePSArecords) {
      out->psarecord.record("id",id);
      out->psarecord.record("state",state);
      out->psarecord.record("ext_grade",ext_grade);
      out->psarecord.record("organised",organised); // only meaningful for mixed_screening
      out->psarecord.record("dx",dx);
      out->psarecord.record("age",age);
      out->psarecord.record("psa",psa);
      out->psarecord.record("t0",t0);
      out->psarecord.record("beta0",beta0);
      out->psarecord.record("beta1",beta1);
      out->psarecord.record("beta2",beta2);
      out->psarecord.record("Z",Z);
    }
    if (!everPSA) {
      if (id < in->nLifeHistories) {
	out->outParameters.revise("age_psa",now());
	// outParameters.revise("first_psa",psa);
      }
      everPSA = true;
    }
    if (formal_costs) {
      add_costs("Invitation");
      lost_productivity(in->panel && psa>=1.0 ? "Formal panel" : "Formal PSA");
      add_costs(in->panel && psa>=1.0 ? "Formal panel" : "Formal PSA");
      scheduleUtilityChange(now(), "Formal PSA");
    } else { // opportunistic costs
      add_costs(in->panel && psa>=1.0 ? "Opportunistic panel" : "Opportunistic PSA");
      lost_productivity(in->panel && psa>=1.0 ? "Opportunistic panel" : "Opportunistic PSA");
      scheduleUtilityChange(now(), "Opportunistic PSA");
    }
    compliance = formal_compliance ?
      in->tableFormalBiopsyCompliance(bounds<double>(psa,3.0,10.0),
				  bounds<double>(age,40,80)) :
      in->tableOpportunisticBiopsyCompliance(bounds<double>(psa,3.0,10.0),
					 bounds<double>(age,40,80));
    // bool positive_test =
    //   (!panel && msg->kind == toScreen && psa >= parameter["psaThreshold"]) ? true :
    //   ( panel && msg->kind == toScreen && biomarker >= parameter["BPThreshold"]) ? true :
    //   (!panel && msg->kind == toBiopsyFollowUpScreen && psa >= parameter["psaThresholdBiopsyFollowUp"]) ? true :
    //   ( panel && msg->kind == toBiopsyFollowUpScreen && biomarker >= parameter["BPThresholdBiopsyFollowUp"]) ? true :
    //   false;
    bool positive_test =
      (msg->kind == toScreen && psa >= in->parameter["psaThreshold"]) ? true :
      (msg->kind == toBiopsyFollowUpScreen && psa >= in->parameter["psaThresholdBiopsyFollowUp"]) ? true :
      false;
    // Important case: PSA<1 (to check)
    // Reduce false positives wrt Gleason 7+ by 1-rFPF: which BPThreshold?
    if (in->panel && positive_test && psa < 10.) {
      if (int(in->parameter("biomarker_model"))==random_correction) { // base model for the biomarker
	if (R::runif(0.0,1.0) < 1.0 - in->parameter["rFPF"])
	  positive_test = false;
      }
      else if (int(in->parameter("biomarker_model"))==psa_informed_correction) { // optimistic model for the biomarker
	if ((ext_grade == ext::Gleason_le_6 &&
	     onset() && psa < in->parameter["PSA_FP_threshold_GG6"]) // FP GG 6 PSA threshold
	    ||  (!onset() && psa < in->parameter["PSA_FP_threshold_nCa"])) {// FP no cancer PSA threshold
	  positive_test = false; // strong assumption
	}
      }
      else {
	REprintf("Parameter biomarker_model not matched: %i\n", int(in->parameter("biomarker_model")));
      }
    }
    if (in->includePSArecords && !onset() && positive_test) {
      out->falsePositives.record("id",id);
      out->falsePositives.record("psa",psa);
      out->falsePositives.record("age",now());
      out->falsePositives.record("age0",t0+35.0);
      out->falsePositives.record("ext_grade",ext_grade);
    }
    // if (panel && !positive_test && t0<now()-35.0 && ext_grade > ext::Gleason_le_6) {
    //   if (R::runif(0.0,1.0) < 1.0-parameter["rTPF"]) positive_test = true;
    // }
    if (positive_test && R::runif(0.0,1.0) < compliance) {
      scheduleAt(now(), toScreenInitiatedBiopsy); // immediate biopsy
    } // assumes similar biopsy compliance, reasonable? An option to different psa-thresholds would be to use different biopsyCompliance. /AK
    else { // re-screening schedules
      if (R::runif(0.0,1.0) < in->parameter["rescreeningCompliance"]) {
	switch (in->screen) {
	case mixed_screening:
	case stockholm3_goteborg:
	case goteborg:
	  {
	    bool mixed_opportunistic = false;
	    if (now() >= in->parameter["start_screening"] && now() < in->parameter["stop_screening"]) {
	      if (psa<1.0 && now()+4.0 <= in->parameter["stop_screening"])
		scheduleAt(now() + 4.0, toScreen);
	      else if (psa>=1.0 && now()+2.0 <= in->parameter["stop_screening"])
		scheduleAt(now() + 2.0, toScreen);
	      else if (in->screen == mixed_screening) {
		opportunistic_rescreening(psa); // older men should start opportunistic rescreening
		mixed_opportunistic = true;
	      }
	      // else do nothing
	    }
	    else if (in->screen == mixed_screening && !mixed_opportunistic) {
	      opportunistic_rescreening(psa); // start opportunistic rescreening
	    }
	  }
	  break;
	case stockholm3_risk_stratified:
	case risk_stratified:
	  if (now() >= in->parameter["start_screening"]) {
	    if (psa<1.0 && now()+8.0 <= in->parameter["stop_screening"])
	      scheduleAt(now() + 8.0, toScreen);
	    if (psa>=1.0 && now()+4.0 <= in->parameter["stop_screening"])
	      scheduleAt(now() + 4.0, toScreen);
	  }
	  break;
	case regular_screen:
	  if (in->parameter["start_screening"] <= now() &&
	      now() + in->parameter["screening_interval"] <= in->parameter["stop_screening"])
	    scheduleAt(now() + in->parameter["screening_interval"], toScreen);
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
	  REprintf("Screening not matched: %s\n",in->screen);
	  break;
	}
      } // rescreening compliance
      if (in->screen == screenUptake || (in->screen == mixed_screening && !organised))
	opportunistic_rescreening(psa); // includes rescreening compliance
    } // rescreening
    in->rngNh->set();
  } break;

  case toClinicalDiagnosis:
    dx = ClinicalDiagnosis;
    RemoveKind(toMetastatic); // competing events
    RemoveKind(toT3plus);
    RemoveKind(toScreen);
    scheduleAt(now(), toClinicalDiagnosticBiopsy); // assumes only one biopsy per clinical diagnosis
    scheduleAt(now(), toTreatment);
    break;

  case toScreenDiagnosis:
    dx = ScreenDiagnosis;
    RemoveKind(toMetastatic); // competing events
    RemoveKind(toT3plus);
    RemoveKind(toClinicalDiagnosis);
    RemoveKind(toScreen);
    scheduleAt(now(), toTreatment);
    break;

  // assumes that biopsies are 100% accurate

  // record additional biopsies for clinical diagnoses
  case toClinicalDiagnosticBiopsy:
    add_costs("Biopsy");
    lost_productivity("Biopsy");
    scheduleUtilityChange(now(), "Biopsy");
    break;

  case toScreenInitiatedBiopsy:
    in->rngScreen->set();
    add_costs("Biopsy");
    lost_productivity("Biopsy");
    scheduleUtilityChange(now(), "Biopsy");

    if (state == Metastatic || (state == Localised && R::runif(0.0, 1.0) < in->parameter["biopsySensitivity"])) { // diagnosed
      scheduleAt(now(), toScreenDiagnosis);
    } else if (!previousNegativeBiopsy) {
      previousNegativeBiopsy = true;
      // first re-screen after negative biopsy
      if (R::runif(0.0,1.0) < in->parameter["rescreeningCompliance"])
	scheduleAt(now() + 1, toBiopsyFollowUpScreen); // schedule one quick PSA retest
    } else if (R::runif(0.0,1.0) < in->parameter["rescreeningCompliance"]) { // next rescreen after negative biopsy
      opportunistic_rescreening(psa); // schedule a routine future screen
    }
    in->rngNh->set();
    break;

  case toTreatment: {
    in->rngTreatment->set();
    double u_tx = R::runif(0.0,1.0);
    double u_adt = R::runif(0.0,1.0);
    if (state == Metastatic) {
      lost_productivity("Metastatic cancer");
      scheduleAt(now(), new cMessageUtilityChange(- in->utility_estimates["Metastatic cancer"]));
    }
    else { // Loco-regional
      tx = calculate_treatment(u_tx,now(),year);
      if (tx == CM) scheduleAt(now(), toCM);
      if (tx == RP) scheduleAt(now(), toRP);
      if (tx == RT) scheduleAt(now(), toRT);
      // check for ADT
      double pADT =
	in->pradt(tx,
	      bounds<double>(now(),50,79),
	      bounds<double>(year,1973,2004),
	      grade);
      if (u_adt < pADT)  {
	adt = true;
	scheduleAt(now(), toADT);
      }
      if (in->debug) Rprintf("id=%i, adt=%d, u=%8.6f, pADT=%8.6f\n",id,adt,u_adt,pADT);
    }
    // reset the random number stream
    in->rngNh->set();
    // check for cure
    bool cured = false;
    double age_c = (state == Localised) ? tc + 35.0 : tmc + 35.0;
    double lead_time = age_c - now();
    // calculate the age at cancer death by c_benefit_type
    double age_cancer_death=R_PosInf;
    if (in->parameter["c_benefit_type"]==LeadTimeBased) { // [new paper ref]
      double pcure = pow(1 - exp(-lead_time * in->parameter["c_benefit_value1"]),
      			 calculate_mortality_hr(age_c));
      if (in->debug) Rprintf("hr for lead time=%f\n", calculate_mortality_hr(age_c));
      cured = (R::runif(0.0,1.0) < pcure);
      if (cured) {
	RemoveKind(toMetastatic);
	RemoveKind(toT3plus);
      } else {
      double u_surv = R::runif(0.0,1.0);
      age_cancer_death = calculate_survival(u_surv,age_c,age_c,calculate_treatment(u_tx,age_c,year+lead_time));
      }
    }
    else if (in->parameter["c_benefit_type"]==StageShiftBased) { // [annals paper ref]
      // calculate survival
      double u_surv = R::runif(0.0,1.0);
      double age_cd = calculate_survival(u_surv,age_c,age_c,calculate_treatment(u_tx,age_c,year+lead_time));
      double age_sd = calculate_survival(u_surv,now(),age_c,tx);
      double weight = exp(- in->parameter["c_benefit_value0"]*lead_time);
      age_cancer_death = weight*age_cd + (1.0-weight)*age_sd;
    }
    else REprintf("c_benefit_type not matched.");
    if (!cured) {
      scheduleAt(age_cancer_death, toCancerDeath);
      // Disutilities prior to a cancer death
      double age_palliative = age_cancer_death - in->utility_duration["Palliative therapy"] - in->utility_duration["Terminal illness"];
      double age_terminal = age_cancer_death - in->utility_duration["Terminal illness"];
      // Reset utilities for those in with a Metatatic diagnosis
      if (state == Metastatic) {
	if (age_palliative > now())
	  scheduleUtilityChange(age_palliative, "Metastatic cancer",false);
	else
	  scheduleUtilityChange(now(), "Metastatic cancer", false);
      }
      if (age_palliative>now()) { // cancer death more than 36 months after diagnosis
	scheduleUtilityChange(age_palliative, "Palliative therapy");
	scheduleUtilityChange(age_terminal, "Terminal illness");
      }
      else if (age_terminal>now()) { // cancer death between 36 and 6 months of diagnosis
	scheduleUtilityChange(now(), "Palliative therapy",false, -1.0);
	scheduleUtilityChange(age_terminal, "Palliative therapy", false, 1.0); // reset
	scheduleUtilityChange(age_terminal,"Terminal illness");
      }
      else // cancer death within 6 months of diagnosis/treatment
	scheduleUtilityChange(now(), "Terminal illness");
    }
    if (in->includeDiagnoses) {
      out->diagnoses.record("id",id);
      out->diagnoses.record("age",age);
      out->diagnoses.record("year",year);
      out->diagnoses.record("psa",psa);
      out->diagnoses.record("ext_grade",ext_grade);
      out->diagnoses.record("ext_state",ext_state);
      out->diagnoses.record("organised",organised); // only meaningful for mixed_screening, keep this?
      out->diagnoses.record("dx",dx);
      out->diagnoses.record("tx",tx);
      out->diagnoses.record("cancer_death",(aoc>age_cancer_death) ? 1.0 : 0.0);
      out->diagnoses.record("age_at_death", (aoc>age_cancer_death) ? age_cancer_death : aoc);
    }
  } break;

  case toRP:
    add_costs("Prostatectomy");
    lost_productivity("Prostatectomy");
    // Scheduling utilities for the first 2 months after procedure
    scheduleUtilityChange(now(), "Prostatectomy part 1");
    // Scheduling utilities for the first 3-12 months after procedure
    scheduleUtilityChange(now() + in->utility_duration["Prostatectomy part 1"],
			  "Prostatectomy part 2");
    break;

  case toRT:
    add_costs("Radiation therapy");
    lost_productivity("Radiation therapy");
    // Scheduling utilities for the first 2 months after procedure
    scheduleUtilityChange(now(), "Radiation therapy part 1");
    // Scheduling utilities for the first 3-12 months after procedure
    scheduleUtilityChange(now() + in->utility_duration["Radiation therapy part 1"],
			  "Radiation therapy part 2");
    break;

  case toCM:
    add_costs("Active surveillance"); // expand here
    lost_productivity("Active surveillance");
    scheduleUtilityChange(now(), "Active surveillance");
    break;

  case toADT:
    // costs & utilities??
    break;

  case toBaselineUtility:
    {
      const cMessageBaselineUtility * msgUtility;
      if ((msgUtility = dynamic_cast<const cMessageBaselineUtility *>(msg)) != 0) {
	baseline_utility = msgUtility->utility;
      } else {
	REprintf("Could not cast to cMessageBaselineUtility.");
      }
    } break;

  case toUtilityChange:
    {
      const cMessageUtilityChange * msgUtilityChange;
      if ((msgUtilityChange = dynamic_cast<const cMessageUtilityChange *>(msg)) != 0) {
	delta_utility += msgUtilityChange->change;
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
  SimInput in;
  SimOutput out;
  
  in.rngNh = new Rng();
  in.rngOther = new Rng();
  in.rngScreen = new Rng();
  in.rngTreatment = new Rng();
  in.rngNh->set();

  // read in the parameters
  List parms(parmsIn);
  List tables = parms["tables"];
  in.parameter = parms["parameter"];
  in.bparameter = parms["bparameter"]; // scalar bools
  List otherParameters = parms["otherParameters"];
  in.debug = as<bool>(parms["debug"]);
  if (! in.bparameter["revised_natural_history"]) {
    in.mubeta2 = as<NumericVector>(otherParameters["mubeta2"]);
    in.sebeta2 = as<NumericVector>(otherParameters["sebeta2"]);
  } else {
    in.mubeta2 = as<NumericVector>(otherParameters["rev_mubeta2"]);
    in.sebeta2 = as<NumericVector>(otherParameters["rev_sebeta2"]);
  }
  NumericVector mu0 = as<NumericVector>(otherParameters["mu0"]);
  in.cost_parameters = as<NumericVector>(otherParameters["cost_parameters"]);
  in.utility_estimates = as<NumericVector>(otherParameters["utility_estimates"]);
  in.utility_duration = as<NumericVector>(otherParameters["utility_duration"]);

  in.production = Table<double,double>(as<DataFrame>(otherParameters["production"]), "ages", "values");
  in.lost_production_proportions = as<NumericVector>(otherParameters["lost_production_proportions"]);

  int n = as<int>(parms["n"]);
  in.includePSArecords = as<bool>(parms["includePSArecords"]);
  in.includeDiagnoses = as<bool>(parms["includeDiagnoses"]);
  int firstId = as<int>(parms["firstId"]);
  in.interp_prob_grade7 =
    NumericInterpolate(as<DataFrame>(tables["prob_grade7"]));
  in.prtxCM = TablePrtx(as<DataFrame>(tables["prtx"]),
			       "Age","DxY","G","CM"); // NB: Grade is now {0,1[,2]} coded cf {1,2[,3]}
  in.prtxRP = TablePrtx(as<DataFrame>(tables["prtx"]),
			       "Age","DxY","G","RP");
  in.pradt = TablePradt(as<DataFrame>(tables["pradt"]),"Tx","Age","DxY","Grade","ADT");
  in.hr_locoregional = TableLocoHR(as<DataFrame>(otherParameters["hr_locoregional"]),"age","ext_grade","psa10","hr");
  in.hr_metastatic = TableMetastaticHR(as<DataFrame>(otherParameters["hr_metastatic"]),"age","hr");
  in.tableOpportunisticBiopsyCompliance = TableBiopsyCompliance(as<DataFrame>(tables["biopsyOpportunisticComplianceTable"]),
						"psa","age","compliance");
  in.tableFormalBiopsyCompliance = TableBiopsyCompliance(as<DataFrame>(tables["biopsyFormalComplianceTable"]),
						"psa","age","compliance");
  in.rescreen_shape = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "shape");
  in.rescreen_scale = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "scale");
  in.rescreen_cure  = TableDDD(as<DataFrame>(tables["rescreening"]), "age5", "total", "cure");

  in.H_dist.clear();
  DataFrame df_survival_dist = as<DataFrame>(tables["survival_dist"]); // Grade,Time,Survival
  DataFrame df_survival_local = as<DataFrame>(tables["survival_local"]); // Age,Grade,Time,Survival
  // extract the columns from the survival_dist data-frame
  IntegerVector sd_grades = df_survival_dist["Grade"];
  NumericVector
    sd_times = df_survival_dist["Time"],
    sd_survivals = df_survival_dist["Survival"];
  for (int i=0; i<sd_grades.size(); ++i)
    in.H_dist[sd_grades[i]].push_back(Double(sd_times[i],-log(sd_survivals[i])));
  for (H_dist_t::iterator it_sd = in.H_dist.begin(); it_sd != in.H_dist.end(); it_sd++)
    it_sd->second.prepare();
  // now we can use: H_dist[grade].invert(-log(u))
  in.H_local.clear();
  // extract the columns from the data-frame
  IntegerVector sl_grades = df_survival_local["Grade"];
  NumericVector
    sl_ages = df_survival_local["Age"],
    sl_times = df_survival_local["Time"],
    sl_survivals = df_survival_local["Survival"];
  // push to the map values and set of ages
  for (int i=0; i<sl_grades.size(); ++i) {
    in.H_local_age_set.insert(sl_ages[i]);
    in.H_local[H_local_t::key_type(sl_ages[i],sl_grades[i])].push_back
      (Double(sl_times[i],-log(sl_survivals[i])));
  }
  // prepare the map values for lookup
  for (H_local_t::iterator it_sl = in.H_local.begin();
       it_sl != in.H_local.end();
       it_sl++)
    it_sl->second.prepare();
  // now we can use: H_local[H_local_t::key_type(*H_local_age_set.lower_bound(age),grade)].invert(-log(u))

  if (in.debug) {
    Rprintf("SurvTime: %f\n",exp(-in.H_local[H_local_t::key_type(65.0,0)].approx(63.934032)));
    Rprintf("SurvTime: %f\n",in.H_local[H_local_t::key_type(* in.H_local_age_set.lower_bound(65.0),0)].invert(-log(0.5)));
    Rprintf("SurvTime: %f\n",exp(-in.H_dist[0].approx(5.140980)));
    Rprintf("SurvTime: %f\n",in.H_dist[0].invert(-log(0.5)));
    // Rprintf("Biopsy compliance: %f\n",tableBiopsyCompliance(bounds<double>(1.0,4.0,10.0), bounds<double>(100.0,55,75)));
    Rprintf("Interp for grade 6/7 (expecting approx 0.3): %f\n",in.interp_prob_grade7.approx(0.143));
    Rprintf("prtxCM(80,2008,1) [expecting 0.970711]: %f\n",in.prtxCM(80.0,2008.0,1));
    {
      double age_diag=51.0;
      ext::grade_t ext_grade = ext::Gleason_ge_8;
      // FhrcPerson person = FhcrcPerson(0,1960);
      // Rprintf("hr_localregional(50,0,)=%g\n",hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext_grade, person.psamean(age_diag)>10 ? 1 : 0));
      Rprintf("hr_localregional(50,8+,0)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext_grade, 0));
      Rprintf("hr_localregional(50,8+,1)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext_grade, 1));
      Rprintf("hr_localregional(50,7,0)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_7, 0));
      Rprintf("hr_localregional(50,7,1)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_7, 1));
      Rprintf("hr_localregional(50,<=6,0)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_le_6, 0));
      Rprintf("hr_localregional(50,<=6,1)=%g\n",in.hr_locoregional(age_diag<50.0 ? 50.0 : age_diag, ext::Gleason_le_6, 1));
      Rprintf("screeningCompliance=%g\n",as<double>(in.parameter["screeningCompliance"]));
    }
  }

  in.nLifeHistories = as<int>(otherParameters["nLifeHistories"]);
  in.screen = as<int>(otherParameters["screen"]);
  if (in.debug) Rprintf("screen=%i\n",in.screen);
  in.panel = as<bool>(parms["panel"]);
  NumericVector cohort = as<NumericVector>(parms["cohort"]); // at present, this is the only chuck-specific data


  // set up the parameters
  double ages0[106];
  boost::algorithm::iota(ages0, ages0+106, 0.0);
  in.rmu0 = Rpexp(&mu0[0], ages0, 106);
  vector<double> ages(101);
  boost::algorithm::iota(ages.begin(), ages.end(), 0.0);
  ages.push_back(1.0e+6);

  // re-set the output objects
  out.report.clear();
  out.shortReport.clear();
  out.costs.clear();
  out.outParameters.clear();
  out.lifeHistories.clear();
  out.psarecord.clear();
  out.falsePositives.clear();
  out.diagnoses.clear();

  out.report.discountRate = in.parameter["discountRate.effectiveness"];
  out.report.setPartition(ages);
  out.shortReport.discountRate = in.parameter["discountRate.effectiveness"];
  out.shortReport.setPartition(ages);
  out.costs.discountRate = in.parameter["discountRate.costs"];
  out.costs.setPartition(ages);

  // main loop
  FhcrcPerson person(&in, &out, 1, 2000);
  for (int i = 0; i < n; ++i) {
    person = FhcrcPerson(&in, &out, i+firstId, cohort[i]);
    Sim::create_process(&person);
    Sim::run_simulation();
    Sim::clear();
    in.rngNh->nextSubstream();
    in.rngOther->nextSubstream();
    in.rngScreen->nextSubstream();
    in.rngTreatment->nextSubstream();
    R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
  }

  // output
  // TODO: clean up these objects in C++ (cf. R)
  return List::create(_("costs") = out.costs.wrap(),                // CostReport
		      _("summary") = out.report.wrap(),             // EventReport
		      _("shortSummary") = out.shortReport.wrap(),   // EventReport
		      _("lifeHistories") = wrap(out.lifeHistories), // vector<LifeHistory::Type>
		      _("parameters") = out.outParameters.wrap(),   // SimpleReport<double>
		      _("psarecord")=out.psarecord.wrap(),          // SimpleReport<double>
		      _("falsePositives")=out.falsePositives.wrap(),// SimpleReport<double>
		      _("diagnoses")=out.diagnoses.wrap()           // SimpleReport<double>
		      );
}

} // anonymous namespace
