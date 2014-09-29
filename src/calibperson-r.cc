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

 Hypothetical microsimulation model. Edited by Alexandra Jauhiainen 130205.

 TODO
 * other causes of death - incorporate rates
 * other transitions
 * age-specific reporting of state probabilities
*/

/**#include "event-r.h"
#include <Rcpp.h>
#include <vector>*/
#include "microsimulation.h"
#include <Rcpp.h>

namespace {

using namespace ssim;
using namespace std;

//! enum of type of disease stage
enum stage_t {DiseaseFree,Precursor,PreClinical,Clinical,Death};

//! enum of type of event type
enum event_t {toPrecursor, toPreClinical, toClinical, toDeath, Count};

//! names of the stages
string stage_names[5] = {"DiseaseFree","Precursor","PreClinical","Clinical","Death"};

//! declare the random number generator
  Rng * rng;

//! Class to simulate a person
class CalibPerson : public cProcess 
{
public:
  stage_t stage;
  bool diseasepot;
  double Lam1,sigm1,p2,lam2,mu3,tau3;
  double clinTime;
  int id;
	
  // static member(s)
  static std::map<std::string, std::vector<double> > report;
  
  static void resetPopulation ();

  CalibPerson() {} // default constructor
  
  CalibPerson(double *par, int i=0) {
    Lam1=par[0];
    sigm1=par[1];
    p2=par[2];
    lam2=par[3];
    mu3=par[4];
    tau3=par[5];
    id=i;
  };
	
  void init();
  virtual void handleMessage(const cMessage* msg);
  virtual Time age() { return now(); }
};

void CalibPerson::resetPopulation() {
  report.clear();
}

// initialise static member(s)
std::map<std::string, std::vector<double> > CalibPerson::report;

/** 
 Initialise a simulation run for an individual
*/
void CalibPerson::init() {
  if (R::runif(0,1)<p2){
    diseasepot = true;}
  else {
    diseasepot = false;
  }
  clinTime=1000;
  stage=DiseaseFree;
  double lam1 = exp(R::rnorm(Lam1,sigm1));
  scheduleAt(R::rexp(lam1), toPrecursor);
  double x = R::runif(0,1);
  scheduleAt((65 - 15*log(-log(x))), toDeath); //Gumbel
  for(int i=10;i<110;i=i+10){
    scheduleAt(i,Count);
  }
}	

/** 
    Handle receiving self-messages
 */
void CalibPerson::handleMessage(const cMessage* msg) {

  double ctime[] = {20,40,60,80};
  int cind;
	
	
  if (msg->kind == toDeath) {
    stage=Death;
    clinTime=std::min(clinTime,now());
	
    for(unsigned int i=0; i<4 ; i++){
      if(i < report["TimeAtRisk"].size()){
	report["TimeAtRisk"][i] += std::min(ctime[i],clinTime);
      }
      else {
	report["TimeAtRisk"].push_back(std::min(ctime[i],clinTime));
      }
      
      if(clinTime < ctime[i]){
	break;
      }
      
    }
	  
    Sim::stop_simulation();
  }
 	
  else if (msg->kind == toPrecursor) {
    stage = Precursor;
    if (diseasepot){
      simtime_t dwellTime = now()+ R::rexp(lam2);
      scheduleAt(dwellTime, toPreClinical);
    }	  
  }
  
  else if (msg->kind == toPreClinical) {
    stage=PreClinical;
    simtime_t dwellTime = now()+ exp(R::rnorm(mu3,tau3*mu3)); 
    scheduleAt(dwellTime, toClinical);
  }
  	
  else if (msg->kind == toClinical) {
    stage=Clinical;
    clinTime = now();
    string stagestr = stage_names[stage];
  }
	   
  else if (msg->kind == Count){
    cind = min(9,int(now()/10 - 1));
    string stagestr = stage_names[stage];
    
    if(report.find(stagestr) == report.end()){ //key not found
      report[stagestr].assign(10,0);
     }
    report[stagestr][cind]+=1;
  }		   
}


extern "C" {

  RcppExport SEXP callCalibrationSimulation(SEXP parms) {
    Rcpp::List parmsl(parms);
    int nin = Rcpp::as<int>(parmsl["n"]);
    std::vector<double> par = Rcpp::as<std::vector<double> >(parmsl["runpar"]);
	  
    CalibPerson::resetPopulation();
    rng = new Rng();
    rng->set();
    
    CalibPerson::report.insert(make_pair("TimeAtRisk", std::vector<double>()));
    
    CalibPerson person;
    for (int i = 0; i < nin; i++) {
      person = CalibPerson(&par[0],0);
      rng->nextSubstream();
      Sim::create_process(&person);
      Sim::run_simulation();
      Sim::clear();
    }

    delete rng;
    
    return Rcpp::wrap(CalibPerson::report);
    
  } 
  
} 

}
