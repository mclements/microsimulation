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
 SSIM and the OMNET++ API. This is specialised for use as an R package
 (#includes and REprintf).

 It also provides several utility classes: Means for statistical
 collection and Rpexp for piecewise constant exponential random number
 generation. It also provides a utility function rweibullHR().

*/

#ifndef MICROSIMULATION_H
#define MICROSIMULATION_H

#include <RcppCommon.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <vector>
#include <utility>
#include <map>

namespace Rcpp {

  // vectors of pairs: recursively defined
  template <class T1, class T2>
    SEXP wrap(const std::vector<std::pair<T1,T2> > v);

  // vectors tuples: enumerated cases
  template <class T1, class T2>
    SEXP wrap(const std::vector<boost::tuple<T1,T2> > v);
  
  template <class T1, class T2, class T3>
    SEXP wrap(const std::vector<boost::tuple<T1,T2,T3> > v);
  
  template <class T1, class T2, class T3, class T4>
    SEXP wrap(const std::vector<boost::tuple<T1,T2,T3,T4> > v);
  
  template <class T1, class T2, class T3, class T4, class T5>
    SEXP wrap(const std::vector<boost::tuple<T1,T2,T3,T4,T5> > v);

  // maps defined in terms of vectors
  template <class T1, class T2>
    SEXP wrap_map(const std::map<T1,T2> v);

}

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Random.h>

#include <siena/ssim.h>
#include "RngStream.h"
#include "rcpp_table.h"

#include <string>
#include <algorithm>
#include <map>
#include <functional>
#include <set>

#include <boost/bind.hpp>
#include <boost/functional.hpp>

namespace ssim {

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::greater;

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
  short kind;
  string name;
  Time sendingTime, timestamp;
  // this does NOT include schedulePriority
  cMessage(const short k = -1, const string n = "") : kind(k), name(n) {
    sendingTime = Sim::clock();
  }
  // currently no setters (keep it lightweight?)
  short getKind() { return kind; }
  string getName() { return name; }
  Time getTimestamp() { return timestamp; }
  Time getSendingTime() {return sendingTime; }
};

/**
   @brief cProcess class for OMNET++ API compatibility.
   This provides a default for Process::process_event() that calls
   cProcess::handleMessage(). This class also provides scheduleAt()
   methods for insert cMessages into the process event queue.
 */
class cProcess : public ssim::Process {
public:
 cProcess() : previousEventTime(0.0) { }; 
  virtual void handleMessage(const cMessage * msg) = 0;
  virtual void process_event(const ssim::Event * e) { // virtual or not?
    const cMessage * msg;
    // cf. static_cast?
    if ((msg = dynamic_cast<const cMessage *>(e)) != 0) { 
      handleMessage(msg);
      previousEventTime = Sim::clock();
    } else {
      // cf. cerr, specialised for R
      REprintf("cProcess is only written to receive cMessage events\n");
    }
  }
  virtual void scheduleAt(Time t, cMessage * msg) { // virtual or not?
    msg->timestamp = t;
    Sim::self_signal_event(msg, t - Sim::clock());  
  }
  virtual void scheduleAt(Time t, string n) {
    scheduleAt(t, new cMessage(-1,n));  
  }
  virtual void scheduleAt(Time t, short k) {
    scheduleAt(t, new cMessage(k));  
  }
  Time previousEventTime;
};

inline bool cMessagePred(const ssim::Event* e, boost::function<bool(const cMessage * msg)> pred) {
    const cMessage * msg = dynamic_cast<const cMessage *>(e);
    return (msg != 0 && pred(msg));
  }

inline bool cMessageNamePred(const ssim::Event* e, const string s) {
    const cMessage * msg = dynamic_cast<const cMessage *>(e);
    return (msg != 0 && msg->name == s); 
  }

inline bool cMessageKindPred(const ssim::Event* e, const short k) {
    const cMessage * msg = dynamic_cast<const cMessage *>(e);
    return (msg != 0 && msg->kind == k); 
  }

/**
   @brief RemoveKind is a function to remove messages with the given kind from the queue (NB: void)
*/
inline void RemoveKind(short kind) {
  return Sim::remove_event(boost::bind(cMessageKindPred,_1,kind));
}

/**
   @brief RemoveName is a function to remove messages with the given name from the queue (NB: void)
*/
inline void RemoveName(string name) {
  return Sim::remove_event(boost::bind(cMessageNamePred,_1,name));
}

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
  Rpexp() {} // blank default constructor
  Rpexp(double *hin, double *tin, int nin) : n(nin) {
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
  }
 /* Rpexp(vector<double> hin, vector<double> tin) : h(hin), t(tin) { */
 /*    n = h.size(); */
 /*    H.resize(n); */
 /*    H[0] = 0.0; */
 /*    if (n>1) { */
 /*      for(i=1;i<n;i++) { */
 /* 	H[i] = H[i-1]+(t[i]-t[i-1])*h[i-1]; */
 /*      } */
 /*    } */
 /*  } */
  double rand(double u, double from = 0.0) {
    double v = 0.0, H0 = 0.0, tstar = 0.0;
    int i = 0, i0 = 0;
    if (from > 0.0) {
      i0 = (from >= t[n-1]) ? (n-1) : int(lower_bound(t.begin(), t.end(), from) - t.begin())-1;
      H0 = H[i0] + (from - t[i0])*h[i0];
    }
    v = -log(u) + H0;
    i = (v >= H[n-1]) ? (n-1) : int(lower_bound(H.begin(), H.end(), v) - H.begin())-1;
    tstar = t[i]+(v-H[i])/h[i];
    return tstar;
  }

 private:
  vector<double> H, h, t;
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
  ~Rng();
  void seed(const unsigned long seed[6]) {
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
  //else if (start == 0.0) return (1.0 - (1.0+discountRate)^(-end)) / log(1.0+discountRate);
  else return (pow(1.0+discountRate,-start) - pow(1.0+discountRate,-end)) / log(1.0+discountRate);
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
template< class Tstate, class Tevent, class T >
class EventReport {
public:
  typedef std::set<T, std::greater<T> > Partition;
  typedef typename Partition::iterator Iterator;
  void setPartition(const vector<T> v) {
    copy(v.begin(), v.end(), inserter(_partition, _partition.begin()));
  }
  void clear() {
    _pt.clear();
    _events.clear();
    _prev.clear();
    _partition.clear();
  }
  void add(const Tstate state, const Tevent event, const T lhs, const T rhs) {
    Iterator lo, hi, it, last;
    lo = _partition.lower_bound(lhs);
    hi = _partition.lower_bound(rhs); 
    last = _partition.begin(); // because it's ordered by greater<T>
    ++_events[boost::make_tuple(state,event,*hi)];
    bool iterating;
    for(it = lo, iterating = true; iterating; iterating = (it != hi), --it) { // decrement for greater<T>
      if (lhs<=(*it) && (*it)<rhs) // cadlag
    	++_prev[std::make_pair(state,*it)];
      if (it == last) {
	_pt[std::make_pair(state,*it)] += rhs - std::max<T>(lhs,*it);
      }
      else {
	T next_value = *(--it); it++; // decrement/increment for greater<T>
	_pt[std::make_pair(state,*it)] += std::min<T>(rhs,next_value) - std::max<T>(lhs,*it);
      }
    }
  }

  SEXP out() {

    using namespace Rcpp;

    return List::create(_("pt") = wrap_map(_pt),
  			_("events") = wrap_map(_events),
  			_("prev") = wrap_map(_prev));
  }

  Partition _partition;
  map<pair<Tstate,T>, int > _prev;
  map<pair<Tstate,T>, T > _pt;
  map<boost::tuple<Tstate,Tevent,T>, int > _events;

};


// http://en.cppreference.com/w/cpp/algorithm/iota
template<class ForwardIterator, class T>
void myiota(ForwardIterator first, ForwardIterator last, T value)
{
    while(first != last) {
        *first++ = value;
        ++value;
    }
}

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
   SEXP wrap(const vector<boost::tuple<T1,T2> > v) {
  int i, n=v.size();
  vector<T1> v1(n);
  vector<T2> v2(n);
  typename vector<boost::tuple<T1,T2> >::const_iterator it;
  for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
    v1[i] = boost::get<0>(*it);
    v2[i] = boost::get<1>(*it);
  }
  return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2));
}

template <class T1, class T2, class T3>
SEXP wrap(const vector<boost::tuple<T1,T2,T3> > v) {
  int i, n=v.size();
  vector<T1> v1(n);
  vector<T2> v2(n);
  vector<T3> v3(n);
  typename vector<boost::tuple<T1,T2,T3> >::const_iterator it;
  for (it=v.begin(), i=0; it<v.end(); ++it, ++i) { 
    v1[i] = boost::get<0>(*it);
    v2[i] = boost::get<1>(*it);
    v3[i] = boost::get<2>(*it);
  }
  return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),_("Var3")=wrap(v3));
}


template <class T1, class T2, class T3, class T4>
SEXP wrap(const vector<boost::tuple<T1,T2,T3,T4> > v) {
  int i, n=v.size();
  vector<T1> v1(n);
  vector<T2> v2(n);
  vector<T3> v3(n);
  vector<T4> v4(n);
  typename vector<boost::tuple<T1,T2,T3,T4> >::const_iterator it;
  for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = boost::get<0>(*it);
      v2[i] = boost::get<1>(*it);
      v3[i] = boost::get<2>(*it);
      v4[i] = boost::get<3>(*it);
  }
  return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
		      _("Var3")=wrap(v3),_("Var4")=wrap(v4));
  }

template <class T1, class T2, class T3, class T4, class T5>
SEXP wrap(const vector<boost::tuple<T1,T2,T3,T4,T5> > v) {
  int i, n=v.size();
  vector<T1> v1(n);
  vector<T2> v2(n);
  vector<T3> v3(n);
  vector<T4> v4(n);
  vector<T5> v5(n);
  typename vector<boost::tuple<T1,T2,T3,T4,T5> >::const_iterator it;
  for (it=v.begin(), i=0; it<v.end(); ++it, ++i) {
      v1[i] = boost::get<0>(*it);
      v2[i] = boost::get<1>(*it);
      v3[i] = boost::get<2>(*it);
      v4[i] = boost::get<3>(*it);
      v5[i] = boost::get<4>(*it);
  }
  return List::create(_("Var1")=wrap(v1),_("Var2")=wrap(v2),
		      _("Var3")=wrap(v3),_("Var4")=wrap(v4),
		      _("Var5")=wrap(v5));
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
   List out =  wrap(x);
   out.push_back(wrap(y),"Value");
   return out;
   //return List::create(_("Key")=wrap(x),_("Value")=wrap(y));
 }

 }


namespace R {
  /**
     @brief rnorm function constrained to be positive. This uses brute-force re-sampling rather
     than conditioning on the distribution function.
  */
  double rnormPos(double mean, double sd);
}

#endif
