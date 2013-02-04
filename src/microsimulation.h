/**
 * @file microsimulation.cc
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

 It also provides two useful utility classes: Means for statistical
 collection and Rpexp for piecewise constant exponential random number
 generation. Lastly, it provides a utility function rweibullHR().

*/

#ifndef EVENT_H
#define EVENT_H

#include <R.h>
#include <Rmath.h>
#include "ssim.h"
#include "Rng_wrapper.h"
#include <string>
#include <vector>

//using namespace std;
//using namespace ssim;
using std::string;
using std::vector;
using ssim::Time;
using ssim::Sim;

// should we use the ssim namespace?

//#define WithRNG(rng,expr) (rng->set(), expr)

/**
   @brief cMessage class for OMNET++ API compatibility.
   This provides a heavier message class than Sim::Event, with
   short 'kind' and string 'name' attributes. 
   The events by default are scheduled using cProcess::scheduleAt(),
   and handled using cProcess::handleMessage() (as per OMNET++).
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
    if ((msg = dynamic_cast<const cMessage *>(e)) != 0) { // cf. static_cast
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


/** 
    Function struct used to compare a message name with a given string
*/
struct cMessageNameEq : public std::binary_function<const ssim::Event *,string,bool> {
  bool operator()(const ssim::Event* e, const string s) const;
};
inline bool cMessageNameEq::operator() (const ssim::Event * e, const string s) const
{ 
  const cMessage * msg = dynamic_cast<const cMessage *>(e);
  return (msg != NULL && msg->name == s); 
};

/** 
    Function struct used to compare a message kind with a given short
*/
struct cMessageKindEq : public std::binary_function<const ssim::Event *,short,bool> {
  bool operator()(const ssim::Event* e, const short k) const;
};
inline bool cMessageKindEq::operator() (const ssim::Event * e, const short k) const
{ 
  const cMessage * msg = dynamic_cast<const cMessage *>(e);
  return (msg != NULL && msg->kind == k); 
};

/**
   Function to remove messages with the given name from the queue (NB: void)
*/
void remove_name(string name);

/**
   Function to remove messages with the given kind from the queue (NB: void)
*/
void remove_kind(short kind);

/**
   simtime_t typedef for OMNET++ API compatibility
*/
typedef Time simtime_t;

/**
   simTime() function for OMNET++ API compatibility
*/
Time simTime();

/**
   now() function for compatibility with C++SIM
*/
Time now();

/**
   @brief Utility class to incrementally add values to calculate the mean,
   sum, variance and standard deviation.
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
   Random number generator class for piecewise constant hazards
 */
class Rpexp {
public: 
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
  double rand(double from = 0.0) {
    double v = 0.0, H0 = 0.0, tstar = 0.0;
    int i = 0, i0 = 0;
    if (from > 0.0) {
      i0 = (from >= t[n-1]) ? (n-1) : int(lower_bound(t.begin(), t.end(), from) - t.begin())-1;
      H0 = H[i0] + (from - t[i0])*h[i0];
    }
    v = rexp(1.0) + H0;
    i = (v >= H[n-1]) ? (n-1) : int(lower_bound(H.begin(), H.end(), v) - H.begin())-1;
    tstar = t[i]+(v-H[i])/h[i];
    return tstar;
  }

 private:
  vector<double> H, h, t;
  int n;
};


/** 
    Random Weibull distribution for a given shape, scale and hazard ratio
*/
double rweibullHR(double shape, double scale, double hr);

#endif