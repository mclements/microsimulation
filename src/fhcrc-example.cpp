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
  
  enum state_t {Healthy,Localised,Metastatic};

  enum diagnosis_t {NotDiagnosed,ClinicalDiagnosis,ScreenDiagnosis};
  
  enum event_t {toLocalised,toMetastatic,toClinicalDiagnosis,toCancerDeath,toOtherDeath,toScreen, 
		toBiopsy,toScreenDiagnosis,toOrganised,toTreatment,toCM,toRP,toRT,toADT};

  enum screen_t {noScreening, randomScreen50to70, twoYearlyScreen50to70, fourYearlyScreen50to70, 
		 screen50, screen60, screen70, screenUptake, stockholm3_goteborg, stockholm3_risk_stratified};

  enum treatment_t {no_treatment, CM, RP, RT};

  typedef boost::tuple<short,short,short,bool,double> FullState;
  //string astates[] = {"stage", "ext_grade", "dx", "psa_ge_3", "cohort"};
  //vector<string> states(astates,astates+5);
  EventReport<FullState,short,double> report;
  CostReport<string> costs;
  map<string, vector<double> > lifeHistories;  // NB: wrap re-defined to return a list
  map<string, vector<double> > parameters;

  bool debug = false;

  double tau2 = 0.0829,
    g0=0.0005, gm=0.0004, gc=0.0015, 
    thetac=19.1334,
    mubeta0=-1.609, sebeta0=0.2384,
    mubeta1=0.04463, sebeta1=0.0430,
    mubeta2[2]={0.0397,0.1678},sebeta2[2]={0.0913,0.3968};

  double c_txlt_interaction = 1.0,
    c_baseline_specific = 1.0;

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

  // initialise input parameters (see R::callFhcrc for actual defaults)
  double screeningCompliance = 0.75;
  double studyParticipation = 35.0/260.0;
  int nLifeHistories = 10, screen = 0;
  double psaThreshold = 3.0;

  // new parameters (we need to merge the old and new implementations)
  double c_low_grade_slope=-0.006;

  Rng * rngNh, * rngOther, * rngScreen, * rngTreatment;
  Rpexp rmu0;

  /** 
      Utilities to record information in a map<string, vector<double> > object
  */
  void record(map<string, vector<double> > & obj, const string variable, const double value) {
    obj[variable].push_back(value);
  }
  void revise(map<string, vector<double> > & obj, const string variable, const double value) {
    *(--obj[variable].end()) = value;
  }


  template<class T>
  T bounds(T x, T a, T b) {
    return (x<a)?a:((x>b)?b:x);
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
    base::grade_t grade;
    ext::grade_t ext_grade;
    treatment_t tx;
    bool adt; 
    double txhaz;
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
  everPSA = previousNegativeBiopsy = organised = adt = false;
  rngNh->set();
  t0 = sqrt(2*R::rexp(1.0)/g0);
  y0 = ymean(t0);
  grade = (R::runif(0.0, 1.0)>=1+c_low_grade_slope*t0) ? base::Gleason_ge_8 : base::Gleason_le_7;
  beta0 = R::rnorm(mubeta0,sebeta0); 
  beta1 = R::rnormPos(mubeta1,sebeta1); 
  beta2 = R::rnormPos(mubeta2[grade],sebeta2[grade]); 
  tm = (log((beta1+beta2)*R::rexp(1.0)/gm + y0) - beta0 + beta2*t0) / (beta1+beta2);
  ym = ymean(tm);
  tc = (log((beta1+beta2)*R::rexp(1.0)/gc + y0) - beta0 + beta2*t0) / (beta1+beta2);
  tmc = (log((beta1+beta2)*R::rexp(1.0)/(gc*thetac) + ym) - beta0 + beta2*t0) / (beta1+beta2);
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
  if (R::runif(0.0,1.0)<screeningCompliance) {
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
  if (R::runif(0.0,1.0)<studyParticipation &&
      ((screen == stockholm3_goteborg) || (screen == stockholm3_risk_stratified)) && 
      (2013.0-cohort>=50.0 && 2013.0-cohort<70.0)) {
    scheduleAt(R::runif(2013.0,2015.0) - cohort, toOrganised);
  }

  rngNh->set();

  // record some parameters
  // faster: initialise the length of the vectors and use an index
  if (id<nLifeHistories) {
    record(parameters,"id",double(id));
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
    record(parameters,"ext_grade",ext_grade);
    record(parameters,"age_psa",-1.0);
    record(parameters,"pca_death",0.0);
  }
}

  // struct Parameters {
  //   int id;
  //   double beta0, beta1, beta2, t0, tm, tc, tmc, y0, ym, aoc, cohort;
  // };
  // typedef boost::tuple<int,double,double,double,double,double,double,double,double,double,double,double> Parameters;
  // std::vector<Parameters> parameters;
  // // Each time, we would do the following:
  // parameters.push_back(Parameters(id,beta0, beta1, beta2, t0, tm, tc, tmc, y0, ym, aoc, cohort));

  // List wrap_parameter(std::vector<Parameters> parameters) {
  //   int n = parameters.size();
  //   vector<int> ids(n);
  //   vector<vector<double> > others(n);
  //   Parameters row;
  //   vector<double> doubles;
  //   for (int i=0; i<n; i++) {
  //     row = parameters[i];
  //     ids[i] = row.id;
  //     others[i] = row.beta0;
  //     others[i].push_back(beta1); // etc
  //   }
  //   return List(ids,others);
  // }

/** 
    Handle self-messages received
 */
void FhcrcPerson::handleMessage(const cMessage* msg) {

  // declarations
  double psa = y(now()-35.0);
  double age = now();
  double year = now() + cohort;

  // record information
  // if (id<10)
  report.add(FullState(state, ext_grade, dx, psa>=3.0, cohort), msg->kind, previousEventTime, now());

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
    if (id<nLifeHistories) {
      record(parameters,"age_d",now());
      revise(parameters,"pca_death",1.0);
    }
    Sim::stop_simulation();
    break;

  case toOtherDeath: 
    if (id<nLifeHistories) 
      record(parameters,"age_d",now());
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
    if (!everPSA) {
      if (id<nLifeHistories)
	revise(parameters,"age_psa",now());
      everPSA = true;
    } 
    costs.add("PSA",now(),200);
    if (psa>=psaThreshold) {
      scheduleAt(now(), toBiopsy); // immediate biopsy
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
    u = pow(u,1/c_baseline_specific); // global improvement to baseline survival
    double lead_time = (tc+35.0) - now();
    double txbenefit = exp(log(txhaz)+log(c_txlt_interaction)*lead_time);
    if (debug) Rprintf("lead_time=%f, txbenefit=%f, u=%f, ustar=%f\n",lead_time,txbenefit,u,pow(u,1/(txbenefit*sxbenefit)));
    u = pow(u,1/(txbenefit*sxbenefit));
    // TODO: calculate survival
    double age_cancer_death;
    if (state == Localised)
      age_cancer_death = tc + 35.0 + H_local[H_local_t::key_type(*H_local_age_set.lower_bound(bounds<double>(now(),50.0,80.0)),grade)].invert(-log(u));
    if (state == Metastatic)
      age_cancer_death = tmc + 35.0 + H_dist[grade].invert(-log(u));
    scheduleAt(age_cancer_death, toCancerDeath);
    if (debug) Rprintf("SurvivalTime=%f, u=%f\n",age_cancer_death -now(), u);
  } break;

  case toRP:
    break;

  case toRT:
    break;

  case toCM:
    break; 

  case toADT:
    break;

  default:
    REprintf("No valid kind of event: %i\n",msg->kind);
    break;
    
  } // switch

} // handleMessage()


RcppExport SEXP callFhcrc(SEXP parmsIn) {

  // declarations
  FhcrcPerson person;

  ssim::RNGScope scope;
  rngNh = new Rng();
  rngOther = new Rng();
  rngScreen = new Rng();
  rngTreatment = new Rng();
  rngNh->set();

  // read in the parameters
  List parms(parmsIn);
  List tables = parms["tables"];
  int n = as<int>(parms["n"]);
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

  nLifeHistories = as<int>(parms["nLifeHistories"]);
  screen = as<int>(parms["screen"]);
  screeningCompliance = as<double>(parms["screeningCompliance"]);
  studyParticipation = as<double>(parms["studyParticipation"]);
  psaThreshold = as<double>(parms["psaThreshold"]);
  NumericVector cohort = as<NumericVector>(parms["cohort"]); // at present, this is the only chuck-specific data

  // set up the parameters
  double ages0[106];
  boost::algorithm::iota(ages0, ages0+106, 0.0);
  rmu0 = Rpexp(mu0, ages0, 106);
  vector<double> ages(101);
  boost::algorithm::iota(ages.begin(), ages.end(), 0.0);
  ages.push_back(1.0e+6);

  // re-set the output objects
  report.clear();
  costs.clear();
  parameters.clear();
  lifeHistories.clear();

  report.setPartition(ages);
  costs.setPartition(ages);

  // main loop
  for (int i = 0; i < n; i++) {
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
		      _("parameters") = wrap(parameters)
		      );
} 

} // anonymous namespace
