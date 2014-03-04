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
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION 

 Microsimulation of prostate cancer using the FHCRC model.

 TODO
 * other causes of death - incorporate rates
 * other transitions
 * age-specific reporting of state probabilities
*/

#include "microsimulation.h"

namespace person_r {

using namespace std;

//! enum for type of Gleason score
enum gleason_t {nogleason,gleasonLt7,gleason7,gleasonGt7};

//! enum of type of disease stage
enum stage_t {Healthy,Localised,DxLocalised,LocallyAdvanced,DxLocallyAdvanced,
	      Metastatic,DxMetastatic,Death};

//! enum of type of event type
enum event_t {toDeath, toPCDeath, toLocalised, toDxLocalised,
	      toDxLocallyAdvanced,
	      toLocallyAdvanced, toMetastatic, toDxMetastatic};

//! Class to simulate a person
class Person : public cProcess 
{
public:
  gleason_t gleason;
  stage_t stage;
  bool dx;
  int id;
  
  // static member(s)
  static std::map<std::string, std::vector<double> > report;
  static std::map<std::string, Rng *> rng;

  static void resetPopulation ();
  //
  Person(int i = 0) : gleason(nogleason), stage(Healthy), dx(false), id(i) {};
  void init();
  virtual void handleMessage(const cMessage* msg);
  virtual Time age() { return now(); }
};

void Person::resetPopulation() {
  report.clear();
}

// initialise static member(s)
std::map<std::string, std::vector<double> > Person::report;
std::map<std::string, Rng *> Person::rng;

/** Hazard ratio for diagnosis
    @param stage Disease stage
 */
double dxHR(stage_t stage) {
  // raise error if healthy?
  return stage==Healthy ? -1 : 
    (stage==Localised ? 1.1308 : 
     (stage==LocallyAdvanced ? 0.5900 :1.3147));
}

/** Hazard ratio for progression
    @param gleason Gleason category
 */
double progressionHR(gleason_t gleason) {
  return gleason==gleasonLt7 ? 1 :
      (gleason==gleason7 ? 1.3874 : 1.4027 * 1.3874);
}

/** 
    Initialise a simulation run for an individual
 */
void Person::init() {
  rng["NH"]->set();
  if (R::runif(0.0,1.0)<0.2241) 
    scheduleAt(R::rweibull(exp(2.3525),64.0218),toLocalised);
  scheduleAt(R::rexp(80.0),toDeath);
}

/** 
    Handle receiving self-messages
 */
void Person::handleMessage(const cMessage* msg) {

  double pDx;

  report["id"].push_back(id);
  report["startTime"].push_back(previousEventTime);
  report["endTime"].push_back(now());
  report["state"].push_back(stage);
  report["event"].push_back(msg->kind);

  rng["NH"]->set(); 

  if (msg->kind == toDeath) {
    Sim::stop_simulation();
  }
 
  else if (msg->kind == toPCDeath) {
    Sim::stop_simulation(); 
  }
  
  else if (msg->kind == toLocalised) {
    stage = Localised;
    gleason = (R::runif(0.0,1.0)<0.6812) ? gleasonLt7 : 
      ((R::runif(0.0,1.0)<0.5016) ? gleason7 : gleasonGt7);
    Time dwellTime = now()+
      rweibullHR(exp(1.0353),19.8617,progressionHR(gleason)*
	       dxHR(stage));
    // now separate out for different transitions
    pDx = 1.1308/(2.1308);
    if (R::runif(0.0,1.0)<pDx) {
      scheduleAt(dwellTime, toDxLocalised);
    }
    else {
      scheduleAt(dwellTime,toLocallyAdvanced);
    }
  }
  
  else if (msg->kind == toLocallyAdvanced) {
    stage=LocallyAdvanced;
    Time dwellTime = now()+
      rweibullHR(exp(1.4404),16.3863,progressionHR(gleason)*
	       dxHR(stage));
    // now separate out for different transitions
    pDx = 0.5900/(1.0+0.5900);
    if (R::runif(0.0,1.0)<pDx) {
      scheduleAt(dwellTime, toDxLocallyAdvanced);
    }
    else {
      scheduleAt(dwellTime,toMetastatic);
    }
  }

  else if (msg->kind == toMetastatic) {
    stage=Metastatic;
    Time dwellTime = now()+
      rweibullHR(exp(1.4404),1.4242,progressionHR(gleason)*
	       dxHR(stage));
    // now separate out for different transitions
    pDx = 1.3147/(1.0+1.3147); 
    if (R::runif(0.0,1.0)<pDx) {
      scheduleAt(dwellTime, toDxMetastatic);
    }
    else {
      scheduleAt(dwellTime,toPCDeath); // prior to diagnosis!
    }
  }

  else if (msg->kind == toDxLocalised) {
    dx=true;
    // relative survival
  }

  else if (msg->kind == toDxLocallyAdvanced) {
    dx=true;
    // relative survival
  }

  else if (msg->kind == toDxMetastatic) {
    dx=true;
    // relative survival
  };
  
};


extern "C" {

  RcppExport SEXP callPersonSimulation(SEXP inseed, SEXP parms) {
    Rcpp::List parmsl(parms);
    Rcpp::IntegerVector inseed2(inseed);
    int nin = Rcpp::as<int>(parmsl["n"]);
    unsigned long seed[6];
    for (int i=0; i<6; i++) {
      seed[i]=(unsigned long)inseed2[i];
    }
    //r_create_current_stream();
    RngStream_SetPackageSeed(seed);
    Person::resetPopulation();
    Person::rng["NH"] = new Rng();
    Person::rng["S"] = new Rng();
    Person::rng["NH"]->set();
    Person person;
    for (int i = 0; i < nin; i++) {
      //Person::rng.foreach(nextSubStream);
      Person::rng["NH"]->nextSubstream();
      Person::rng["S"]->nextSubstream();
      person = Person(i);
      Sim::create_process(&person);
      Sim::run_simulation();
      Sim::clear();
    }
    // tidy up -- what needs to be deleted?
    delete Person::rng["NH"];
    delete Person::rng["S"];
    // Person::rng.clear();
    // output arguments to R
    return Rcpp::wrap(Person::report);
    
  } // callPersonSimulation()
  
} // extern "C"

} // namespace person_r

namespace {

  class VerySimple : public cProcess {
  public:
    void init() {
      scheduleAt(10.0, "a message");
      scheduleAt(11.0, "another message");
    };
    virtual void handleMessage(const cMessage* msg) {};
  };
  
extern "C" {

  RcppExport SEXP callSpeedTest() {
    VerySimple simple;
    for (int i = 0; i < 1000000; i++) {
      simple = VerySimple();      
      Sim::create_process(&simple);
      Sim::run_simulation();
      Sim::clear();
    }
    return Rcpp::wrap(1);
    
  } // callSpeedTest()
  
} // extern "C"

} // anonymous namespace
    


