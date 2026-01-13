/**
 * @file microsimulation.h
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
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION

 cMessage and cProcess classes, providing some compatability between
 SSIM and the OMNET++ API. This is specialised for use as an R package.

 It also provides several utility classes: Means for statistical
 collection and Rpexp for piecewise constant exponential random number
 generation. It also provides a utility function rweibullHR().

*/

#ifndef MICROSIMULATION_H
#define MICROSIMULATION_H

#include <RcppCommon.h>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <utility>
#include <map>
#include <string>

/* https://stackoverflow.com/questions/7110301/generic-hash-for-tuples-in-unordered-map-unordered-set for a generic hash function for std::tuple
 */
#include <tuple>
namespace std{
  namespace
  {
    // Code from boost
    // Reciprocal of the golden ratio helps spread entropy
    //     and handles duplicates.
    // See Mike Seymour in magic-numbers-in-boosthash-combine:
    //     http://stackoverflow.com/questions/4948780
    template <class T>
    inline void hash_combine(std::size_t& seed, T const& v)
    {
      seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }
    // Recursive template code derived from Matthieu M.
    template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
    struct HashValueImpl
    {
      static void apply(size_t& seed, Tuple const& tuple)
      {
	HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
	hash_combine(seed, std::get<Index>(tuple));
      }
    };
    template <class Tuple>
    struct HashValueImpl<Tuple,0>
    {
      static void apply(size_t& seed, Tuple const& tuple)
      {
	hash_combine(seed, std::get<0>(tuple));
      }
    };
  }
  template <typename ... TT>
  struct hash<std::tuple<TT...>>
  {
    size_t
    operator()(std::tuple<TT...> const& tt) const
    {
      size_t seed = 0;
      HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
      return seed;
    }
  };
  template <typename T1, typename T2>
  struct hash<std::pair<T1,T2>>
  {
    size_t
    operator()(std::pair<T1,T2> const& p) const
    {
      size_t seed = 0;
      hash_combine(seed, p.first);
      hash_combine(seed, p.second);
      return seed;
    }
  };
}

namespace Rcpp {

  // vectors of pairs: recursively defined
  template <class T1, class T2>
    SEXP wrap(const std::vector<std::pair<T1,T2> > v);

  // vectors tuples: enumerated cases
  template <class T1, class T2>
    SEXP wrap(const std::vector<std::tuple<T1,T2> > v);

  template <class T1, class T2, class T3>
    SEXP wrap(const std::vector<std::tuple<T1,T2,T3> > v);

  template <class T1, class T2, class T3, class T4>
    SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4> > v);

  template <class T1, class T2, class T3, class T4, class T5>
    SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4,T5> > v);

  template <class T1, class T2, class T3, class T4, class T5, class T6>
    SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4,T5,T6> > v);

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7>
    SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4,T5,T6,T7> > v);

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8>
    SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8> > v);

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
    SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> > v);

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9, class T10>
    SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > v);

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9, class T10, class T11>
  SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11> > v);

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9, class T10, class T11, class T12>
  SEXP wrap(const std::vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12> > v);

  // maps defined in terms of vectors
  template <class T1, class T2>
    SEXP wrap_map(const std::map<T1,T2> v);

  template <class T1a, class T1b, class T2>
    SEXP wrap_map(const std::map<std::pair<T1a,T1b>,T2> v,
		  std::string key, std::string name1, std::string name2);

  template <class T1a, class T1b, class T1c, class T2>
    SEXP wrap_map(const std::map<std::tuple<T1a,T1b,T1c>,T2> v,
		  std::string key, std::string name1, std::string name2, std::string name3);

  template <class T1, class T2, class T3, class T4, class T5>
  SEXP wrap_map(const std::map<std::pair<std::tuple<T1,T2,T3>,T4>,T5> v,
		  std::string key, std::string name1, std::string name2);

  // unordered_maps defined in terms of vectors
  template <class T1, class T2>
    SEXP wrap_map(const std::unordered_map<T1,T2> v);

  template <class T1a, class T1b, class T2>
    SEXP wrap_map(const std::unordered_map<std::pair<T1a,T1b>,T2> v,
		  std::string key, std::string name1, std::string name2);

  template <class T1a, class T1b, class T1c, class T2>
    SEXP wrap_map(const std::unordered_map<std::tuple<T1a,T1b,T1c>,T2> v,
		  std::string key, std::string name1, std::string name2, std::string name3);

  template <class T1, class T2, class T3, class T4, class T5>
    SEXP wrap_map(const std::unordered_map<std::pair<std::tuple<T1,T2,T3>,T4>,T5> v,
		  std::string key, std::string name1, std::string name2);


} // namespace Rcpp

#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Random.h>

#include <siena/ssim.h>
#include <RngStream.h>
#include <rcpp_table.h>
#include <gsm.h>

#include <string>
#include <algorithm>
#include <functional>
#include <set>

namespace ssim {

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::greater;

// forward declarations
 class cProcess;

/**
   @brief WithRNG is a macro for using the random number generator rng
   and then evaluating the expression expr.
*/
#define WithRNG(rng,expr) (rng->set(), expr)

/**
   @brief cMessage class for OMNET++ API compatibility.  This provides
   a heavier message class than Sim::Event, with short 'kind' and
   std::string 'name' attributes.  The events by default are scheduled
   using cProcess::scheduleAt(), and handled using
   cProcess::handleMessage() (as per OMNET++).  NB:
   cProcess::scheduleAt() uses simulation time rather than time in
   state (which is used by Sim::self_signal_event()).
*/
class cMessage : public ssim::Event {
public:
  cMessage(const short k = -1, const string n = "", const ProcessId process_id = -1,
	   const ProcessId sender_process_id = -1, const short priority = 0) :
    ssim::Event(), kind(k), name(n), sendingTime(-1.0), timestamp(0), process_id(process_id),
    sender_process_id(sender_process_id) {
    this->priority = priority;
  }
  // currently no setters (keep it lightweight?)
  short getKind() { return kind; }
  string getName() { return name; }
  Time getTimestamp() { return timestamp; }
  Time getSendingTime() {return sendingTime; }
  short getSchedulePriority() {return priority; }
  short kind;
  string name;
  Time sendingTime, timestamp;
  ProcessId process_id, sender_process_id;
  string str() const {
    std::ostringstream stringStream;
    stringStream << "kind=";
    stringStream << kind;
    stringStream << ",name=";
    stringStream << name;
    string str = stringStream.str();
    return str;
  }
};

/**
   @brief Cancel is a function to cancel messages who satisfy a predicate.

   This works across all processes.
*/
  inline void
  Cancel(std::function<bool(const cMessage * msg)> pred) {
    return Sim::ignore_event([pred](const Event * e) {
			       const cMessage * msg = dynamic_cast<const cMessage *>(e);
			       return (msg != 0 && pred(msg));
			     });
  }

/**
   @brief RemoveKind is a function to remove messages with the given kind from the queue

   This works across all processes.
*/
  inline void RemoveKind(short kind) {
    Cancel([kind](const cMessage * msg) { return msg->kind == kind; });
  }

/**
   @brief CancelKind is a function to remove messages with the given kind from the queue

   This works across all processes.
*/
  inline void CancelKind(short kind) {
    Cancel([kind](const cMessage * msg) { return msg->kind == kind; });
  }

 /**
    @brief RemoveName is a function to remove messages with the given name from the queue

   This works across all processes.
 */
  inline void RemoveName(string name) {
    Cancel([name](const cMessage * msg) { return msg->name == name; });
  }
 /**
    @brief CancelName is a function to remove messages with the given name from the queue

   This works across all processes.
 */
  inline void CancelName(string name) {
    Cancel([name](const cMessage * msg) { return msg->name == name; });
  }

  /**
     @brief CancelEvents is a function to remove messages from the queue

     This works across all processes.
  */
inline void CancelEvents() {
    Cancel([](const cMessage * msg) { return true; });
  }


/**
   @brief cProcess class for OMNET++ API compatibility.

   This provides implementations for the ProcessWithPId including:

   <ol>
   <li> process_event(const Event*) that calls handleMessage(const cMessage*)
   <li> initialize() that calls init().
   </ol>

   The user can also implement a stop() method which is called at the
   end of a process. The scheduleAt() and send() methods insert
   cMessages into the current and other process event queues,
   respectively. Messages for this process can be cancelled using
   cancel(), cancel_kind(), cancel_name() and cancel_events(). The
   previous event time for this process can be found using previous().

   Events are ordered by scheduled time and priority (NB: smaller priority goes earlier).
   
 */
class cProcess : public ssim::ProcessWithPId {
public:
  Time startTime, previousEventTime;
  cProcess(Time startTime = Time(0.0)) : startTime(startTime), previousEventTime(startTime) { }
  /**
      @brief Abstract method to handle each message
   */
  virtual void handleMessage(const cMessage * msg) = 0;
  /**
      @brief Abstract method to initialise a given process.
   */
  virtual void init() = 0;
  // see also `void stop()`
  /**
      @brief schedules at time t a message msg to the current process.
      Adds the sendingTime, process_id and sender_process_id to the message.
   */
  virtual void scheduleAt(Time t, cMessage * msg, short priority=0) { // virtual or not?
    msg->timestamp = t;
    msg->sendingTime = Sim::clock();
    msg->process_id = msg->sender_process_id = pid();
    msg->priority = priority;
    Sim::self_signal_event(msg, t - Sim::clock());
  }
  /**
      @brief schedules at time t a message msg with a specific name to the current process.
   */
  virtual void scheduleAt(Time t, string name, short priority=0) {
    scheduleAt(t, new cMessage(-1,name), priority);
  }
  /**
      @brief schedules at time t a message msg with a specific kind to the current process.
   */
  virtual void scheduleAt(Time t, short k, short priority=0) {
    scheduleAt(t, new cMessage(k,""), priority);
  }
  /**
      @brief send to a given process at time t a message msg.
      Adds the sendingTime, process_id and sender_process_id to the message.
   */
  virtual void send(ProcessId process_id, Time t, cMessage * msg, short priority=0) { // virtual or not?
    msg->timestamp = t;
    msg->process_id = process_id;
    msg->sender_process_id = pid();
    msg->sendingTime = Sim::clock();
    msg->priority = priority;
    Sim::signal_event(process_id, msg, t - Sim::clock());
  }
  /**
      @brief sends to a given process at time t a message msg with a specific name.
   */
  virtual void send(ProcessId process_id, Time t, string name, short priority=0) {
    send(process_id, t, new cMessage(-1,name), priority);
  }
  /**
      @brief sends to a given process at time t a message msg with a specific kind.
   */
  virtual void send(ProcessId process_id, Time t, short kind, short priority=0) {
    send(process_id, t, new cMessage(kind,""), priority);
  }
  /**
      @brief Given a predicate, remove the messages for this process.
   */
  void cancel(std::function<bool(const cMessage * msg)> pred) {
    Cancel([pred,this](const cMessage * msg) {
	     return pred(msg) &&
	       msg->process_id == this->pid();
	   });
  }
  /**
      @brief Given a kind, remove the messages for this process.
   */
  void cancel_kind(short kind) {
    cancel([kind,this](const cMessage * msg) {
	     return msg->kind == kind &&
	       msg->process_id == this->pid();
	   });
  }
  /**
      @brief Given a name, remove the messages for this process.
   */
  void cancel_name(string name) {
    cancel([name,this](const cMessage * msg) {
	     return msg->name == name &&
	       msg->process_id == this->pid();
	   });
  }
  /**
      @brief Remove all of the messages for this process.
   */
  void cancel_events() {
    cancel([this](const cMessage * msg) {
	     return msg->process_id == this->pid();
	   });
  }
  /**
      @brief Returns the value of the previous event time for this process.
   */
  double previous() {
    return previousEventTime;
  }
  /**
      @brief Implements the method required by ssim::ProcessWithPId.
      Calls init() and sets the previousEventTime to the startTime.
   */
  void initialize() { init(); previousEventTime = startTime; }
  /**
      @brief Implements the method required by ssim::ProcessWithPId
      This requires that the event* can be cast to a cMessage*.
      This calls the handleMessage(const cMessage*) method and updates
      previousEventTime.
   */
  virtual void process_event(const ssim::Event * e) { // virtual or not?
    const cMessage * msg;
    if ((msg = dynamic_cast<const cMessage *>(e)) != 0) {
      handleMessage(msg);
      previousEventTime = Sim::clock();
    } else {
      // cf. cerr, specialised for R
      REprintf("cProcess is only written to receive cMessage events\n");
    }
  }
};

/**
   @brief simtime_t typedef for OMNET++ API compatibility
*/
typedef Time simtime_t;

/**
   @brief simTime() function for OMNET++ API compatibility
*/
Time simTime();

/**
   @brief now() function for compatibility with C++SIM
*/
Time now();

/**
   @brief Utility class to incrementally add values to calculate the mean,
   sum, variance and standard deviation. This could be replaced by boost::accumulator.
 */
class Means {
public:
  double mean() { return _sum/_n; }
  double var() {return (long double)_n/(_n-1)*(_sumsq/_n - mean()*mean()); }
  int n() {return _n;}
  double sum() {return _sum;}
  double sd() {return sqrt(var()); }
  Means() : _n(0), _sum(0.0), _sumsq(0.0) {}
  Means* operator +=(const double value) {
    _n++;
    _sum += (long double) value;
    _sumsq += (long double)value * (long double)value;
    return this;
  }
  SEXP wrap() {
    return Rcpp::DataFrame::create(_("n")=n(),
				   _("mean")=mean(),
				   _("var")=var(),
				   _("sd")=sd(),
				   _("se")=sd()/sqrt(n()),
				   _("sum")=sum(),
				   _("sumsq")=_sumsq
				   );
  }
  //friend std::ostream& operator<<(std::ostream& os, Means& m);
private:
  int _n;
  long double _sum, _sumsq;
};

/**
   @brief Rpexp is a random number generator class for piecewise constant hazards.
   Given time lower bounds t and piecewise constant hazards h, rand() returns a random time.
   The random number is calculated using the inversion formula.
   Constructors provided for arrays.
 */
class Rpexp {
public:
  double Z;
  Rpexp() { Z=1.0; } // blank default constructor
  Rpexp(double *hin, double *tin, int nin) : Z(1.0), n(nin) {
    int i;
    H.resize(n);
    t.resize(n);
    h.resize(n);
    H[0]=0.0; h[0]=hin[0]; t[0]=tin[0];
    if (n>1) {
      for(i=1;i<n;i++) {
	h[i]=hin[i]; t[i]=tin[i];
	H[i] = H[i-1]+(t[i]-t[i-1])*h[i-1];
      }
    }
    hi = h;
    Hi = H;
  }
  Rpexp(SEXP hin_, SEXP tin_) {
    Rcpp::NumericVector hin = Rcpp::as<Rcpp::NumericVector>(hin_);
    Rcpp::NumericVector tin = Rcpp::as<Rcpp::NumericVector>(tin_);
    int i;
    Z = 1.0;
    n = hin.size();
    H.resize(n);
    t.resize(n);
    h.resize(n);
    H[0]=0.0; h[0]=hin[0]; t[0]=tin[0];
    if (n>1) {
      for(i=1;i<n;i++) {
	h[i]=hin[i]; t[i]=tin[i];
	H[i] = H[i-1]+(t[i]-t[i-1])*h[i-1];
      }
    }
    hi = h;
    Hi = H;
  }
  void set_Z(double Z) {
    this->Z = Z;
    for(int i=0;i<n;i++) {
	hi[i] = Z*h[i];
	Hi[i] = Z*H[i];
    }
  }
  double life_expectancy(double s, double Z_error = 1.0) {
    double value = 0.0;
    double survival = 1.0;
    for(int j=floor(s); j<(int)h.size(); j++) {
      double prev_survival = survival;
      double delta = double(j)+1-std::max(s,double(j));
      survival *= exp(-hi[j]*delta*Z_error);
      if (prev_survival==survival) value += survival*delta;
      else value += (prev_survival-survival)/(hi[j]*Z_error);
    }
    return value;
  }
  double n_year_risk(double s, double n=5.0, double Z_error = 1.0) {
    double survival = 1.0;
    for(int j=floor(s); j<floor(s+n); j++) {
      double delta = std::min(double(j)+1,s+n)-std::max(s,double(j));
      survival *= exp(-hi[j]*delta*Z_error);
    }
    if (floor(n+s) < n+s) { // add the last bit:)
      int j = floor(s+n);
      double delta = s+n-floor(s+n);
      survival *= exp(-hi[j]*delta*Z_error);
    }
    return 1.0-survival;
  }
  double rand(double u, double from = 0.0) {
    double v = 0.0, H0 = 0.0, tstar = 0.0;
    int i = 0, i0 = 0;
    if (from > 0.0) {
      i0 = (from >= t[n-1]) ? (n-1) : int(lower_bound(t.begin(), t.end(), from) - t.begin())-1;
      H0 = Hi[i0] + (from - t[i0])*hi[i0];
    }
    v = -log(u) + H0;
    i = (v >= Hi[n-1]) ? (n-1) : int(lower_bound(Hi.begin(), Hi.end(), v) - Hi.begin())-1;
    tstar = t[i]+(v-Hi[i])/hi[i];
    return tstar;
  }
 private:
  vector<double> H, h, t, Hi, hi;
  int n;
};


/**
    @brief Random Weibull distribution for a given shape, scale and hazard ratio
*/
double rweibullHR(double shape, double scale, double hr);


/**
    @brief C++ wrapper class for the RngStream library.
    set() sets the current R random number stream to this stream.
    This is compliant with being a Boost random number generator.
*/
static int counter_id = 0;
class Rng : public RngStream {
 public:
  typedef double result_type;
  result_type operator()() { return RandU01(); }
  result_type min() { return 0.0; }
  result_type max() { return 1.0; }
  Rng() : RngStream() { id = ++counter_id; }
  virtual ~Rng();
  void seed(const double seed[6]) {
    SetSeed(seed);
  }
  void set();
  void nextSubstream() { ResetNextSubstream(); }
  int id;
};


extern "C" { // functions that will be called from R

  /**
      @brief A utility function to create the current_stream.
      Used when initialising the microsimulation package in R.
  */
  void r_create_current_stream();

  /**
      @brief A utility function to remove the current_stream.
      Used when finalising the microsimulation package in R.
  */
  void r_remove_current_stream();

  /**
      @brief A utility function to set the user random seed for the simulation.
  */
  void r_set_user_random_seed(double * seed);

  /**
      @brief A utility function to set the user random seed for the simulation.
  */
  void r_get_user_random_seed(double * seed);

  /* /\**  */
  /*     @brief A utility function to move to the next user random stream for the simulation. */
  /* *\/ */
  /* void r_next_rng_stream(); */

  /**
      @brief A utility function to move to the next user random stream for the simulation.
  */
  void r_next_rng_substream();

  /**
      @brief A utility function to advance the random sub-stream n steps for a specified seed.
  */
  void r_rng_advance_substream(double * seed, int * n);

  /**
      @brief Simple test of the random streams (with a stupid name)
  */
  void test_rstream2(double * x);

} // extern "C"


  /**
      @brief Simple function to calculate the integral between the start and end times
      for (1+kappa)^(-u), where kappa is the discountRate (e.g. 0.03)
  */
inline double discountedInterval(double start, double end, double discountRate) {
  if (discountRate == 0.0) return end - start;
  else return (pow(1.0+discountRate,-start) - pow(1.0+discountRate,-end)) / log(1.0+discountRate);
}

  /**
      @brief Simple function to calculate the integral between the start and end times
      for y*(1+kappa)^(-u), where kappa is the discountRate (e.g. 0.03)
  */
  inline double discountedInterval(double y, double start, double end, double discountRate) {
    if (discountRate == 0.0) return y*(end - start);
  else return y*(pow(1.0+discountRate,-start) - pow(1.0+discountRate,-end)) / log(1.0+discountRate);
}

  /**
     @brief Simple function to calculate 1/(1+discountRate)^(time)
  */
  inline double discountedPoint(double time, double discountRate) {
    return discountRate <= 0.0 ? 1.0 : pow(1.0+discountRate,-time);
  }
  /**
     @brief Simple function to calculate y/(1+discountRate)^(time)
  */
  inline double discountedPoint(double y, double time, double discountRate) {
    return discountRate <= 0.0 ? y : y*pow(1.0+discountRate,-time);
  }


/**
    @brief Function to transpose a vector of vectors.
    This assumes that all inner vectors have the same size and
    allocates space for the complete result in advance.
    From  http://stackoverflow.com/questions/6009782/how-to-pivot-a-vector-of-vectors
*/
 template <class T>
   std::vector<std::vector<T> > transpose(const std::vector<std::vector<T> > data) {
   std::vector<std::vector<T> > result(data[0].size(),
				       std::vector<T>(data.size()));
    for (typename std::vector<T>::size_type i = 0; i < data[0].size(); i++)
      for (typename std::vector<T>::size_type j = 0; j < data.size(); j++) {
	result[i][j] = data[j][i];
        }
    return result;
 }


 /**
   @brief EventReport class for collecting statistics on person-time, prevalence and numbers of events.
 */
 template< class State, class Event, class Time = double, class Utility = double>
   class EventReport {
 public:
 typedef std::set<Time, std::greater<Time> > Partition; // NB: greater<Time> cf. lesser<Time>
 typedef typename Partition::iterator Iterator;
 typedef std::pair<State,Time> Pair;
 typedef std::unordered_map<pair<State,Time>, int > PrevMap;
 typedef std::unordered_map<pair<State,Time>, Utility > UtilityMap;
 typedef std::unordered_map<pair<State,Time>, Time > PtMap;
 typedef std::unordered_map<std::tuple<State,Event,Time>, int > EventsMap;
 typedef vector<Utility> IndividualUtilities;
   EventReport(Utility discountRate = 0.0, bool outputUtilities = true, int size = 1, Time startReportAge = Time(0), bool indiv = false) :
     discountRate(discountRate), outputUtilities(outputUtilities), startReportAge(startReportAge), id(0), indiv(indiv) {
   _vector.resize(size);
   setPartition(startReportAge);
 }
 void resize(int size) {
   _vector.resize(size);
 }
 void setPartition(const vector<Time> v) {
   copy(v.begin(), v.end(), inserter(_partition, _partition.begin()));
  }
 void setPartition(const Time start = 0.0, const Time finish = 100.0,
		   const Time delta = 1.0,
		   const Time maxTime = Time(1.0e100)) {
   _partition.clear();
   for (Time t=start; t<=finish; t+=delta) _partition.insert(t);
   _partition.insert(maxTime);
 }
 void setIndivN(const int n) {
   resize(n);
   indiv = true;
 }
 void setStartReportAge(const Time a) {
   startReportAge = a;
 }
 void clear() {
   _ut.clear();
   _pt.clear();
   _events.clear();
   _prev.clear();
   _partition.clear();
   _vector.clear();
   current = Utility(0);
 }
 void individualReset () {
   if (now() >= startReportAge)
     mean_utilities += double(current);
   if (indiv) {
     if (id>=int(_vector.size()))
       REprintf("Vector too small in EventReport: use resize(int) method");
     else
       _vector[id] = (now() >= startReportAge) ? double(current) : NA_REAL;
   }
   current = Utility(0);
   id++;
 }
 Utility discountedUtilities(Time a, Time b, Utility utility = 1.0) {
   if (discountRate == 0.0) return utility * (b-a);
   else if (a==b) return 0.0;
   else if (discountRate>0.0) {
     double alpha = log(1.0+discountRate);
     return utility/alpha*(exp(-a*alpha) - exp(-b*alpha));
   }
   else {
     REprintf("discountRate less than zero.");
     return 0.0;
   }
 }
 void addBrief(const Time lhs, const Time rhs, const Utility utility = 1) {
   if (rhs >= startReportAge)
     current += discountedUtilities(std::max<Time>(lhs-startReportAge,Time(0)), rhs-startReportAge, utility);
 }
 void add(const State state, const Event event, const Time lhs, const Time rhs, const Utility utility = 1, int index = 0 /* deprecated argument */) {
   addBrief(lhs, rhs, utility);
   Iterator lo, hi, it, last;
   lo = _partition.lower_bound(lhs);
   hi = _partition.lower_bound(rhs);
   last = _partition.begin(); // because it is ordered by greater<Time>
   ++_events[std::make_tuple(state,event,*hi)];
   bool iterating;
   for(it = lo, iterating = true; iterating; iterating = (it != hi), --it) { // decrement for greater<Time>
      if (lhs<=(*it) && (*it)<rhs) // cadlag
    	++_prev[Pair(state,*it)];
      if (it == last) {
	if (outputUtilities) {
	  Utility u = discountedUtilities(std::max<Time>(lhs,*it), rhs, utility);
	  _ut[Pair(state,*it)] += u;
	}
	_pt[Pair(state,*it)] += rhs - std::max<Time>(lhs,*it);
      }
      else {
	Time next_value = *(--it); it++; // decrement/increment for greater<Time>
	if (outputUtilities) {
	  Utility u = discountedUtilities(std::max<Time>(lhs,*it), std::min<Time>(rhs,next_value), utility);
	  _ut[Pair(state,*it)] += u;
	}
	_pt[Pair(state,*it)] += std::min<Time>(rhs,next_value) - std::max<Time>(lhs,*it);
      }
   }
 }
 template<class T>
 void append_map(T & base_map, T & new_map) {
   typename T::iterator it;
   for (it = new_map.begin(); it != new_map.end(); ++it)
     base_map[it->first] += it->second;
 }
 void append(EventReport<State,Event,Time,Utility> & er) {
   append_map<PrevMap>(_prev,er._prev);
   append_map<EventsMap>(_events,er._events);
   append_map<PtMap>(_pt,er._pt);
   append_map<UtilityMap>(_ut,er._ut);
   _vector.insert(_vector.end(), er._vector.begin(), er._vector.end());
 }
 SEXP wrap() {
   using namespace Rcpp;
   if (_events.size() == 0) return List::create();
     else if (outputUtilities)
       return List::create(_("pt") = wrap_map<State,Time,Time>(_pt,"Key","age","pt"),
			   _("ut") = wrap_map<State,Time,Utility>(_ut,"Key","age","utility"),
			   _("events") = wrap_map<State,Event,Time,int>(_events,"Key","event","age","number"),
			   _("prev") = wrap_map<State,Time,int>(_prev,"Key","age","number"));
   else
     return List::create(_("pt") = wrap_map<State,Time,Time>(_pt,"Key","age","pt"),
			 _("events") = wrap_map<State,Event,Time,int>(_events,"Key","event","age","number"),
			 _("prev") = wrap_map<State,Time,int>(_prev,"Key","age","number"));
 }
 SEXP wrap_indiv() {
   return Rcpp::wrap(_vector);
 }
 SEXP wrap_means() {
   return mean_utilities.wrap();
 }
 Utility discountRate, current;
 bool outputUtilities;
 Partition _partition;
 PrevMap _prev;
 UtilityMap _ut;
 PtMap _pt;
 EventsMap _events;
 IndividualUtilities _vector;
 Means mean_utilities;
 Time startReportAge;
 int id;
 bool indiv;
};

 /**
   @brief SummaryReport class for collecting statistics on person-time, prevalence, events and costs.
 */
  template< class State, class Event = short, class Time = double, class Utility = double, class Cost = double>
  class SummaryReport {
  public:
    typedef std::set<Time, std::greater<Time> > Partition; // NB: greater<Time> cf. lesser<Time>
    typedef typename Partition::iterator Iterator;
    typedef std::pair<State,Time> Pair;
    typedef std::unordered_map<pair<State,Time>, int > PrevMap;
    typedef std::unordered_map<pair<State,Time>, Utility > UtilityMap;
    typedef std::unordered_map<pair<State,Time>, Time > PtMap;
    typedef std::unordered_map<std::tuple<State,Event,Time>, int > EventMap;
    typedef std::unordered_map<pair<State,Time>, Cost > CostMap;
    typedef std::vector<Cost> IndividualCosts;
    typedef std::vector<Utility> IndividualUtilities;
    SummaryReport(int n = 1, bool indivp = true, Utility discountRate = 0.0) : n(n), indivp(indivp) {
      resize(indivp ? n : 1);
      setDiscountRate(discountRate);
      setUtility(1.0);
      setCost(0.0);
    }
    void resize(int size) {
      _indivUtilities.resize(size);
      _indivCosts.resize(size);
    }
    void setDiscountRate(Utility discountRate) {
      setUtilityDiscountRate(discountRate);
      setCostDiscountRate(discountRate);
    }
    void setUtilityDiscountRate(Utility discountRate) {
      this->utilityDiscountRate = discountRate;
      utilityAlpha = log(1.0+utilityDiscountRate);
    }
    void setCostDiscountRate(Cost discountRate) {
      this->costDiscountRate = discountRate;
      costAlpha = log(1.0+costDiscountRate);
    }
    void setPartition(const vector<Time> v) {
      copy(v.begin(), v.end(), inserter(_partition, _partition.begin()));
    }
    void setPartition(const Time start, const Time finish, const Time delta,
		      const Time maxTime = Time(1.0e100)) {
      _partition.clear();
      for (Time t=start; t<=finish; t+=delta) _partition.insert(t);
      _partition.insert(maxTime);
    }
    void setUtility(Utility _utility) {
      this->utility = _utility;
    }
    void setCost(Cost _cost) {
      this->cost = _cost;
    }
    void clear() {
      _ut.clear();
      _pt.clear();
      _events.clear();
      _prev.clear();
      _partition.clear();
      _costs.clear();
      _indivUtilities.clear();
      _indivCosts.clear();
    }
    void add(const State state, const Event event, const Time lhs, const Time rhs, int index = 0) {
      // Time lhs = process->previousEventTime;
      // Time rhs = ssim::now();
      if (!indivp) index = 0;
      Iterator lo, hi, it, last;
      lo = _partition.lower_bound(lhs);
      hi = _partition.lower_bound(rhs);
      last = _partition.begin(); // because it is ordered by greater<Time>
      ++_events[std::make_tuple(state,event,*hi)];
      bool iterating;
      for(it = lo, iterating = true; iterating; iterating = (it != hi), --it) { // decrement for greater<Time>
	if (lhs<=(*it) && (*it)<rhs) // cadlag
	  ++_prev[Pair(state,*it)];
	if (it == last) {
	  Utility u = discountedUtilityInterval(std::max<Time>(lhs,*it), rhs, utility);
	  _ut[Pair(state,*it)] += u;
	  _indivUtilities[index] += u;
	  Cost c = discountedCostInterval(std::max<Time>(lhs,*it), rhs, cost);
	  _costs[Pair(state,*it)] += c;
	  _indivCosts[index] += c;
	  _pt[Pair(state,*it)] += rhs - std::max<Time>(lhs,*it);
	}
	else {
	  Time next_value = *(--it); it++; // decrement/increment for greater<Time>
	  Utility u = discountedUtilityInterval(std::max<Time>(lhs,*it), std::min<Time>(rhs,next_value), utility);
	  _ut[Pair(state,*it)] += u;
	  _indivUtilities[index] += u;
	  Cost c = discountedCostInterval(std::max<Time>(lhs,*it), std::min<Time>(rhs,next_value), cost);
	  _costs[Pair(state,*it)] += c;
	  _indivCosts[index] += c;
	  _pt[Pair(state,*it)] += std::min<Time>(rhs,next_value) - std::max<Time>(lhs,*it);
	}
      }
    }
    // void addPointCost(const State state, const Time time, const Cost cost, const int index = 0) {
    void addPointCost(const State state, const Cost cost, int index = 0) {
      if (!indivp) index = 0;
      Time time = ssim::now();
      Time time_lhs = * _partition.lower_bound(time);
      Cost c = discountedCost(time,cost);
      _costs[Pair(state,time_lhs)] += c;
      _indivCosts[index] += c;
    }
    void append(SummaryReport<State,Event,Time,Utility> & er) {
      append_map<PrevMap>(_prev,er._prev);
      append_map<EventMap>(_events,er._events);
      append_map<PtMap>(_pt,er._pt);
      append_map<UtilityMap>(_ut,er._ut);
      append_map<CostMap>(_costs,er._costs);
      _indivUtilities.insert(_indivUtilities.end(), er._indivUtilities.begin(), er._indivUtilities.end());
      _indivCosts.insert(_indivCosts.end(), er._indivCosts.begin(), er._indivCosts.end());
    }
    Rcpp::List asList() {
      using namespace Rcpp;
      if (_events.size() == 0) return List::create();
      Rcpp::DataFrame indiv =
	Rcpp::DataFrame::create(_("utilities")=_indivUtilities,
				_("costs")=_indivCosts);
      Rcpp::List out = List::create(_("pt") = wrap_map<State,Time,Time>(_pt,"Key","age","pt"),
				    _("ut") = wrap_map<State,Time,Utility>(_ut,"Key","age","utility"),
				    _("events") = wrap_map<State,Event,Time,int>(_events,"Key","event","age","number"),
				    _("prev") = wrap_map<State,Time,int>(_prev,"Key","age","number"),
				    _("costs") = wrap_map<State,Time,Cost>(_costs,"Key","age","cost"),
				    _("n")=n,
				    _("indivp")=indivp,
				    _("utilityDiscountRate")=utilityDiscountRate,
				    _("costDiscountRate")=costDiscountRate,
				    _("indiv")=indiv);
      out.attr("class")="SummaryReport";
      return out;
    }
    Utility discountedUtilityInterval(Time a, Time b, Utility utility) {
      if (a == b || utility == 0.0) return 0.0;
      else if (utilityDiscountRate == 0.0) return utility * (b-a);
      else if (utilityDiscountRate>0.0) {
	return utility/utilityAlpha*(exp(-a*utilityAlpha) - exp(-b*utilityAlpha));
      }
      else {
	REprintf("utilityDiscountRate less than zero: %f.\n", utilityDiscountRate);
	return 0.0;
      }
    }
    Cost discountedCostInterval(Time a, Time b, Cost cost) {
      if (a == b || cost == 0.0) return 0.0;
      else if (costDiscountRate == 0.0) return cost * (b-a);
      else if (costDiscountRate>0.0) {
	return cost/costAlpha*(exp(-a*costAlpha) - exp(-b*costAlpha));
      }
      else {
	REprintf("costDiscountRate less than zero: %f.\n", costDiscountRate);
	return 0.0;
      }
    }
    Cost discountedCost(Time a, Cost cost) {
      if (cost == 0.0) return 0.0;
      else if (costDiscountRate == 0.0) return cost;
      else if (costDiscountRate>0)
	return cost*exp(-a*costAlpha);
      else {
	REprintf("costDiscountRate less than zero: %f.\n", costDiscountRate);
	return 0;
      }
    }
    template<class T>
    void append_map(T & base_map, T & new_map) {
      typename T::iterator it;
      for (it = new_map.begin(); it != new_map.end(); ++it)
	base_map[it->first] += it->second;
    }
    int n;
    bool indivp;
    Partition _partition;
    PrevMap _prev;
    UtilityMap _ut;
    PtMap _pt;
    EventMap _events;
    CostMap _costs;
    IndividualUtilities _indivUtilities;
    IndividualCosts _indivCosts;
    // cProcess *process;
    Utility utilityDiscountRate, utilityAlpha, utility;
    Cost costDiscountRate, costAlpha, cost;
  };

 /**
    @brief CostReport class for collecting statistics on costs.
 */
 template< class State, class Time = double, class Cost = double>
   class CostReport {
 public:
 typedef std::set<Time, std::greater<Time> > Partition;
 typedef std::pair<State,Time> Pair;
 typedef CostReport<State,Time,Cost> This;
 typedef std::unordered_map<pair<State,Time>, Cost > Table;
 typedef std::vector<Cost> IndividualCosts;
   CostReport(Cost discountRate = 0, int size = 1, Time startReportAge = Time(0), bool indiv = false) : discountRate(discountRate), startReportAge(startReportAge), id(0), indiv(indiv) {
   _vector.resize(size);
 }
 void individualReset () {
   if (now() >= startReportAge)
     mean_costs += double(current);
   if (indiv) {
     if (id>=int(_vector.size()))
       REprintf("Vector too small in CostReport: use resize(int) method");
     else
       _vector[id] = (now() >= startReportAge) ? double(current) : NA_REAL;
   }
   current = Cost(0);
   id++;
 }
 void setIndivN(const int n) {
   resize(n);
   indiv = true;
 }
 void setStartReportAge(const Time a) {
   startReportAge = a;
 }
 Cost discountedCost(Time a, Cost cost) {
   if (discountRate == 0.0) return cost;
   else if (discountRate>0.0)
     return cost/pow(1+discountRate,a);
   else {
     REprintf("discountRate less than zero.");
     return 0;
   }
 }
 void setPartition(const vector<Time> v) {
   copy(v.begin(), v.end(), inserter(_partition, _partition.begin()));
 }
 void setPartition(const Time start, const Time finish, const Time delta,
		   const Time maxTime = Time(1.0e100)) {
   _partition.clear();
   for (Time t=start; t<=finish; t+=delta) _partition.insert(t);
   _partition.insert(maxTime);
 }
 void clear() {
   _table.clear();
   _partition.clear();
   _vector.clear();
   current = Cost(0);
 }
 void resize(int size) {
   _vector.resize(size);
 }
 void append(This & new_report) { // assuming that discountRate and _partition are the same for both reports
   typename Table::iterator it;
   _vector.insert(_vector.start(), new_report._vector.begin(), new_report._vector.end());
   for(it = new_report._table.begin(); it != new_report._table.end(); ++it) {
     _table[it->first] += it->second;
   }
 }
   // // This class is currently only suitable for point costs:(
   // void setCost(const State state, const Cost cost) {
   //   _cost = cost;
   // }
   void addPointCost(const State state, const Time time, const Cost cost, const int index = 0) {
     add(state, time, cost, index);
   }
   // note that this is equivalent to SummaryReport::addPointCost() with state information
   void add(const State state, const Time time, const Cost cost, const int index = 0 /* deprecated argument */) {
   Time time_lhs = * _partition.lower_bound(time);
   Cost c = discountedCost(time,cost);
   _table[Pair(state,time_lhs)] += c;
   if (startReportAge == 0.0)
     current += c;
   else
     if (time>=startReportAge)
       current += discountedCost(time-startReportAge,cost);
 }
 SEXP wrap() {
   return Rcpp::wrap_map(_table,"Key","age","cost");
 }
 SEXP wrap_indiv() {
   return Rcpp::wrap(_vector);
 }
 SEXP wrap_means() {
   return mean_costs.wrap();
 }
 Cost discountRate, current;
 Partition _partition;
 Table _table;
 IndividualCosts _vector;
 Means mean_costs;
 Time startReportAge;
 int id;
 bool indiv;
 };

 /**
    @brief SimpleReport class for collecting data for homogeneous fields of type T with string names.
 */
 template<class T = double>
   class SimpleReport {
 public:
 /**
    @brief Map typedef for a map of strings to vector<T>
 */
 typedef map<string,vector<T> > Map;
 /**
    @brief record adds a value for a given field into the SimpleReport
 */
 void record(string field, T value) {
   _data[field].push_back(value);
 }
 /**
    @brief revise changes the last value in a given field.
    Reminder: std::map::operator[] creates an element for a key if it does not exist.
 */
 void revise(string field, T value) {
   if (!_data[field].empty())
     _data[field].pop_back();
   _data[field].push_back(value);
 }
 /**
    @brief clear the report data
 */
 void clear() { _data.clear(); }
 /**
    @brief wrap the report as a DataFrame or a List
 */
 SEXP wrap() {
   return Rcpp::wrap(_data);
 }
 /**
    @brief append another SimpleReport, which is useful for aggregating multiple reports.
 */
 void append(SimpleReport<T> & obj) {
   for(typename Map::iterator it = obj._data.begin(); it != obj._data.end(); ++it) {
     _data[it->first].insert(_data[it->first].end(), it->second.begin(), it->second.end());
   }
 }
 /**
    @brief _data class member of a map from strings to vector<T>.
 */
 Map _data;
 };

  double R_zeroin2(			/* An estimate of the root */
					double ax,				/* Left border | of the range	*/
					double bx,				/* Right border| the root is seeked*/
					double fa, double fb,		/* f(a), f(b) */
					double (*f)(double x, void *info),	/* Function under investigation	*/
					void *info,				/* Add'l info passed on to f	*/
					double *Tol,			/* Acceptable tolerance		*/
					int *Maxit);				/* Max # of iterations */

  /** 
      Adapt a function object (functor) to work with R_zeroin2()
  **/
  template<class Functor>
  double R_zeroin2_adaptor(double x, void * par) {
    Functor * functor = (Functor *) par;
    return functor->operator()(x);
  }

  /** 
      Use R_zeroin2 with a function object (functor)
      @return tuple<double,double,int> with (root, Tol, Maxit)
  **/
  template<class Functor>
  std::tuple<double,double,int>
  R_zeroin2_functor(double a, double b, Functor functor, double eps = 1.0e-8) {
    double Tol = eps;
    int Maxit = 100;
    double root = R_zeroin2(a,b,functor(a),functor(b),&R_zeroin2_adaptor<Functor>,(void *) &functor,
			    &Tol, &Maxit);
    return std::make_tuple(root, Tol, Maxit);
  }

  /** 
      Use R_zeroin2 with a function object pointer (functor)
      @return tuple<double,double,int> with (root, Tol, Maxit)
  **/
  template<class Functor>
  std::tuple<double,double,int>
  R_zeroin2_functor_ptr(double a, double b, Functor *functor, double tol = 1.0e-8, int maxit = 100) {
    double Tol = tol;
    int Maxit = maxit;
    double root = R_zeroin2(a,b,(*functor)(a),(*functor)(b),
			    &R_zeroin2_adaptor<Functor>,(void *) functor,
			    &Tol, &Maxit);
    return std::make_tuple(root, Tol, Maxit);
  }

  // MVK distribution
  // Reminder:
  // nu: initiation rate
  // alpha: cell division rate
  // beta: cell death rate
  // mu: malignant transformation rate
  // A=(beta+mu-alpha - sqrt((beta+mu-alpha)^2+4*mu*alpha))/2
  // B=(beta+mu-alpha + sqrt((beta+mu-alpha)^2+4*mu*alpha))/2
  // g=-(A+B)=alpha-beta-mu is approximately equal to the net proliferation rate
  // B: upper bound for the malignant transformation rate
  // delta=nu/alpha
  // B-A = sqrt((beta+mu-alpha)^2+4*mu*alpha)
  //
  // Bounds:
  // (nu, alpha, beta, mu) >= 0
  // delta>0, B>=0, A<=0 => B-A>=B>=0
  inline double dMVK(double t, double A, double B, double delta) {
    double P = expm1((B-A) * t) * exp(B * delta * t) * pow(B - A, delta);
    double Q = pow(B * exp((B - A) * t) - A, 1+delta);
    double val = - delta * A * B * P/Q;
    if (!R_finite(val)) val = 0.0;
    return val;
  }
  inline double pMVK(double t, double A, double B, double delta, bool lower_tail = true) {
    double logS = delta*(log(B-A) + B*t - log(B*exp((B-A)*t) - A));
    return lower_tail ? -expm1(logS) : exp(logS);
  }
  inline double qMVK_fast_unstable(double u, double A, double B, double delta,
		     double x0 = 100.0, int iter = 20) {
    for (int i=0; i<iter; i++) {
      x0 = x0 - (pMVK(x0, A, B, delta) - u) / dMVK(x0, A, B, delta);
    }
    return x0;
  }
  inline double qMVK(double u, double A, double B, double delta,
		     double upper = 1e4) {
    int maxiter = 10;
    double ubound = upper, out = 0.0;
    auto functor = [&](double x) { return pMVK(x, A, B, delta) - u; };
    for (int i = 0; i<maxiter; i++) { 
      std::tuple<double,double,int> values =
	R_zeroin2_functor(0.0, ubound, functor, 1.0e-8);
      out = std::get<0>(values);
      if (out>0 && out<ubound) break;
      if (i<maxiter-1) ubound *= 10.0;
      }
    if (out==0 || out==ubound) out = NA_REAL;
    return out;
  }
  inline double rMVK(double u, double A, double B, double delta,
		     double upper = 1e4) {
    return qMVK(u, A, B, delta, upper);
  }
  inline double hMVK(double t, double A, double B, double delta) {
    return dMVK(t,A,B,delta)/pMVK(t,A,B,delta,false);
  }

  /** 
      Class for a multivariate normal distribution.
  **/
  class MVN {
  public:
    arma::vec mu;
    arma::mat Sigma, L;
    MVN() {}; // null default constructor -- I hope that's okay
    MVN(arma::vec mu, arma::mat Sigma) : mu(mu), Sigma(Sigma) {
      L = arma::chol(Sigma, "lower");
    }
    arma::vec rand() {
      arma::vec z(mu.size());
      for (size_t j=0; j<mu.size(); j++) z[j] = R::rnorm(0.0,1.0);
      return mu + arma::vec(L * z);
    }
    arma::mat randn(size_t n = 1) {
      arma::mat y(n,mu.size());
      for (size_t i=0; i<n; i++) {
	arma::vec u = rand();
	for (size_t j=0; j<mu.size(); j++) y(i,j) = u[j];
      }
      return y;
    }
  };
  // //[[Rcpp::export]]
  // arma::mat test_MVN2(size_t n, arma::vec mu, arma::mat Sigma) {
  //   return ssim::MVN(mu,Sigma).randn(n);
  // }
  
} // namespace ssim

namespace Rcpp {

  template <class T1, class T2>
    SEXP wrap(const vector<pair<T1,T2> > v) {
    vector<T1> v1;
    vector<T2> v2;
    typename vector<pair<T1,T2> >::const_iterator it;
    for (it=v.begin(); it<v.end(); ++it) {
      v1.push_back(it->first);
      v2.push_back(it->second);
    }
   return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2));
  }

  template <class T1, class T2>
    SEXP wrap(const vector<std::tuple<T1,T2> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    typename vector<std::tuple<T1,T2> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2));
  }

  template <class T1, class T2, class T3>
    SEXP wrap(const vector<std::tuple<T1,T2,T3> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    typename vector<std::tuple<T1,T2,T3> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),_("Var3")=wrap(v3));
  }


  template <class T1, class T2, class T3, class T4>
    SEXP wrap(const vector<std::tuple<T1,T2,T3,T4> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    typename vector<std::tuple<T1,T2,T3,T4> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4));
  }

  template <class T1, class T2, class T3, class T4, class T5>
    SEXP wrap(const vector<std::tuple<T1,T2,T3,T4,T5> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    vector<T5> v5(n);
    typename vector<std::tuple<T1,T2,T3,T4,T5> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
      v5[i] = std::get<4>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4),
			_("Var5")=wrap(v5));
  }

  template <class T1, class T2, class T3, class T4, class T5, class T6>
    SEXP wrap(const vector<std::tuple<T1,T2,T3,T4,T5,T6> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    vector<T5> v5(n);
    vector<T6> v6(n);
    typename vector<std::tuple<T1,T2,T3,T4,T5,T6> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
      v5[i] = std::get<4>(*it);
      v6[i] = std::get<5>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4),
			_("Var5")=wrap(v5),_("Var6")=wrap(v6));
  }

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7>
    SEXP wrap(const vector<std::tuple<T1,T2,T3,T4,T5,T6,T7> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    vector<T5> v5(n);
    vector<T6> v6(n);
    vector<T7> v7(n);
    typename vector<std::tuple<T1,T2,T3,T4,T5,T6,T7> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
      v5[i] = std::get<4>(*it);
      v6[i] = std::get<5>(*it);
      v7[i] = std::get<6>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4),
			_("Var5")=wrap(v5),_("Var6")=wrap(v6),
			_("Var7")=wrap(v7));
  }

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8>
    SEXP wrap(const vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    vector<T5> v5(n);
    vector<T6> v6(n);
    vector<T7> v7(n);
    vector<T8> v8(n);
    typename vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
      v5[i] = std::get<4>(*it);
      v6[i] = std::get<5>(*it);
      v7[i] = std::get<6>(*it);
      v8[i] = std::get<7>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4),
			_("Var5")=wrap(v5),_("Var6")=wrap(v6),
			_("Var7")=wrap(v7),_("Var8")=wrap(v8));
  }

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
    SEXP wrap(const vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    vector<T5> v5(n);
    vector<T6> v6(n);
    vector<T7> v7(n);
    vector<T8> v8(n);
    vector<T9> v9(n);
    typename vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
      v5[i] = std::get<4>(*it);
      v6[i] = std::get<5>(*it);
      v7[i] = std::get<6>(*it);
      v8[i] = std::get<7>(*it);
      v9[i] = std::get<8>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4),
			_("Var5")=wrap(v5),_("Var6")=wrap(v6),
			_("Var7")=wrap(v7),_("Var8")=wrap(v8),
			_("Var9")=wrap(v9));
  }

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9, class T10>
    SEXP wrap(const vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    vector<T5> v5(n);
    vector<T6> v6(n);
    vector<T7> v7(n);
    vector<T8> v8(n);
    vector<T9> v9(n);
    vector<T10> v10(n);
    typename vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
      v5[i] = std::get<4>(*it);
      v6[i] = std::get<5>(*it);
      v7[i] = std::get<6>(*it);
      v8[i] = std::get<7>(*it);
      v9[i] = std::get<8>(*it);
      v10[i] = std::get<9>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4),
			_("Var5")=wrap(v5),_("Var6")=wrap(v6),
			_("Var7")=wrap(v7),_("Var8")=wrap(v8),
			_("Var9")=wrap(v9),_("Var10")=wrap(v10));
  }
  
  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9, class T10, class T11>
  SEXP wrap(const vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    vector<T5> v5(n);
    vector<T6> v6(n);
    vector<T7> v7(n);
    vector<T8> v8(n);
    vector<T9> v9(n);
    vector<T10> v10(n);
    vector<T11> v11(n);
    typename vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
      v5[i] = std::get<4>(*it);
      v6[i] = std::get<5>(*it);
      v7[i] = std::get<6>(*it);
      v8[i] = std::get<7>(*it);
      v9[i] = std::get<8>(*it);
      v10[i] = std::get<9>(*it);
      v11[i] = std::get<10>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4),
			_("Var5")=wrap(v5),_("Var6")=wrap(v6),
			_("Var7")=wrap(v7),_("Var8")=wrap(v8),
			_("Var9")=wrap(v9),_("Var10")=wrap(v10),
			_("Var11")=wrap(v11));
  }

  template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9, class T10, class T11, class T12>
  SEXP wrap(const vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12> > v) {
    int i, n=v.size();
    vector<T1> v1(n);
    vector<T2> v2(n);
    vector<T3> v3(n);
    vector<T4> v4(n);
    vector<T5> v5(n);
    vector<T6> v6(n);
    vector<T7> v7(n);
    vector<T8> v8(n);
    vector<T9> v9(n);
    vector<T10> v10(n);
    vector<T11> v11(n);
    vector<T12> v12(n);
    typename vector<std::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12> >::const_iterator it;
    for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = std::get<0>(*it);
      v2[i] = std::get<1>(*it);
      v3[i] = std::get<2>(*it);
      v4[i] = std::get<3>(*it);
      v5[i] = std::get<4>(*it);
      v6[i] = std::get<5>(*it);
      v7[i] = std::get<6>(*it);
      v8[i] = std::get<7>(*it);
      v9[i] = std::get<8>(*it);
      v10[i] = std::get<9>(*it);
      v11[i] = std::get<10>(*it);
      v12[i] = std::get<11>(*it);
    }
    return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
			_("Var3")=wrap(v3),_("Var4")=wrap(v4),
			_("Var5")=wrap(v5),_("Var6")=wrap(v6),
			_("Var7")=wrap(v7),_("Var8")=wrap(v8),
			_("Var9")=wrap(v9),_("Var10")=wrap(v10),
			_("Var11")=wrap(v11),_("Var12")=wrap(v12));
  }

  template <class T1, class T2>
    SEXP wrap_map(const std::map<T1,T2> v) {
    int i;
    int n = v.size();
    vector<T1> x(n);
    vector<T2> y(n);
    typename std::map<T1,T2>::const_iterator it;
    for (it=v.begin(), i=0; it != v.end(); ++it, ++i) {
      x[i] = (*it).first;
      y[i] = (*it).second;
    }
    return DataFrame::create(_("Key")=wrap(x),_("Value")=wrap(y));
  }

  // Special cases for wrap_map:
  //
  // map<pair<Tstate,T>, int > _prev;
  // map<pair<Tstate,T>, T > _pt;
  // map<std::tuple<Tstate,Tevent,T>, int > _events;

  template <class T1a, class T1b, class T1c, class T2>
    SEXP wrap_map(const std::map<std::tuple<T1a,T1b,T1c>,T2> v,
		  std::string key, std::string name1, std::string name2, std::string name3) {
    typedef std::tuple<T1a,T1b,T1c> Tuple;
    int i;
    int n = v.size();
    vector<T1a> xa(n);
    vector<T1b> xb(n);
    vector<T1c> xc(n);
    vector<T2> y(n);
    typename std::map<Tuple,T2>::const_iterator it;
    for (it=v.begin(), i=0; it != v.end(); ++it, ++i) {
      xa[i] = get<0>((*it).first);
      xb[i] = get<1>((*it).first);
      xc[i] = get<2>((*it).first);
      y[i] = (*it).second;
    }
    return Rcpp::DataFrame::create(_[key]=Rcpp::wrap(xa), _[name1]=Rcpp::wrap(xb),
				   _[name2]=Rcpp::wrap(xc), _[name3]=Rcpp::wrap(y));
  }


  template <class T1a, class T1b, class T2>
    SEXP wrap_map(const std::map<std::pair<T1a,T1b>,T2> v,
		  std::string key, std::string name1, std::string name2) {
    typedef std::pair<T1a,T1b> Pair;
    int i;
    int n = v.size();
    vector<T1a> xa(n);
    vector<T1b> xb(n);
    vector<T2> y(n);
    typename std::map<Pair,T2>::const_iterator it;
    for (it=v.begin(), i=0; it != v.end(); ++it, ++i) {
      xa[i] = (*it).first.first;
      xb[i] = (*it).first.second;
      y[i] = (*it).second;
    }
    return Rcpp::DataFrame::create(_[key]=Rcpp::wrap(xa), _[name1]=Rcpp::wrap(xb),
				   _[name2]=Rcpp::wrap(y));
   }


  template <class T1, class T2>
    SEXP wrap_map(const std::unordered_map<T1,T2> ov) {
    std::map<T1,T2> v(ov.begin(), ov.end());
    return wrap_map<T1,T2>(v);
  }

  template <class T1a, class T1b, class T1c, class T2>
    SEXP wrap_map(const std::unordered_map<std::tuple<T1a,T1b,T1c>,T2> ov,
		  std::string key, std::string name1, std::string name2, std::string name3) {
    std::map<std::tuple<T1a,T1b,T1c>,T2> v(ov.begin(), ov.end());
    return wrap_map<T1a, T1b, T1c, T2>(v, key, name1, name2, name3);
  }

  template <class T1a, class T1b, class T2>
    SEXP wrap_map(const std::unordered_map<std::pair<T1a,T1b>,T2> ov,
		  std::string key, std::string name1, std::string name2) {
    std::map<std::pair<T1a,T1b>,T2> v(ov.begin(), ov.end());
    return wrap_map<T1a, T1b, T2>(v, key, name1, name2);
   }

  template <class T1, class T2, class T3, class T4, class T5>
  SEXP wrap_map(const std::map<std::pair<std::tuple<T1,T2,T3>,T4>,T5> v,
		std::string name1, std::string name2) {
    typedef std::tuple<T1,T2,T3> Tuple;
    typedef std::pair<Tuple,T4> Pair;
    int i;
    int n = v.size();
    vector<T1> x1(n);
    vector<T2> x2(n);
    vector<T3> x3(n);
    vector<T4> x4(n);
    vector<T5> y(n);
    typename std::map<Pair,T5>::const_iterator it;
    for (it=v.begin(), i=0; it != v.end(); ++it, ++i) {
      Tuple tuple = (*it).first.first;
      x1[i] = std::get<0>(tuple);
      x2[i] = std::get<1>(tuple);
      x3[i] = std::get<2>(tuple);
      x4[i] = (*it).first.second;
      y[i] = (*it).second;
    }
    DataFrame out = Rcpp::DataFrame::create(Rcpp::Named("col1")=wrap(x1),
					    Rcpp::Named("col2")=wrap(x2),
					    Rcpp::Named("col3")=wrap(x3),
					    Rcpp::Named(name1)=wrap(x4),
					    Rcpp::Named(name2)=wrap(y));
    return out;
   }

  template <class T1, class T2, class T3, class T4, class T5>
  SEXP wrap_map(const std::unordered_map<std::pair<std::tuple<T1,T2,T3>,T4>,T5> ov,
		std::string name1, std::string name2) {
    std::map<std::pair<std::tuple<T1,T2,T3>,T4>,T5> v(ov.begin(), ov.end());
    return wrap_map<T1,T2,T3,T4,T5>(v, name1, name2);
   }



} // Rcpp namespace


namespace R {
  /**
     @brief rnorm function constrained to be positive. This uses
     brute-force re-sampling rather than conditioning on the
     distribution function.
  */
  double rnormPos(double mean, double sd);

  /**
     @brief rllogis function for a random covariate from a
     log-logistic distribution with shape and scale.  S(t) =
     1/(1+(t/scale)^shape).
  */
  double rllogis(double shape, double scale);
  /**
     @brief rllogis_trunc function for a random covariate from a
     log-logistic distribution with shape and scale with minimum time
     left.  S(t|t>x)=S(t)/S(x) where S(t)=1/(1+(t/scale)^shape).
  */
  double rllogis_trunc(double shape, double scale, double left);
  /**
     @brief Random draw from a Gompertz distribution.

     Importantly, the parameterisation is as per flexsurv, where
     hazard(t)=rate*exp(shape*t).  As a reminder, if shape<0 then
     there is a non-zero probability of cure (that is, no event),
     which is represented by R_PosInf.
  */
  double rgompertz(double shape = 1.0, double rate = 1.0);
}

#endif
