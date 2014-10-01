#include "microsimulation.h"
#include <omp.h>

namespace {

  using namespace std;
  using namespace Rcpp;
  using namespace ssim;

  enum Event {toOtherDeath};

  class Person : public cProcess 
  {
  public:
    int id;
    Person(const int i = 0) : id(i) { };
    void init();
    virtual void handleMessage(const cMessage* msg); 
  };

  /** 
      Initialise a simulation run for an individual
  */
  void Person::init() {
    double aoc = 85.0; //rmu0.rand(runif(genNh));
    scheduleAt(aoc,toOtherDeath);
    //Rprintf("id=%i started.\n",id);
  }
  
  /** 
      Handle self-messages received
  */
  void Person::handleMessage(const cMessage* msg) {
    
    // handle messages by kind
    
    switch(msg->kind) {
      
    case toOtherDeath:
      sim->stop_simulation();
      //Rprintf("id=%i stopped.\n",id);
      break;
      
    default:
      //REprintf("No valid kind of event: %i\n",msg->kind);
      break;
      
    } // switch
    
  } // handleMessage()
  
  
  RcppExport SEXP callSimTest(SEXP parmsIn) {
    
    List parms(parmsIn);
    int n = as<int>(parms["n"]);
    
    int nthreads, tid, i;

#pragma omp parallel shared(nthreads) private(i,tid)
    {
      /* This section only prints the number of available threads and the start each thread*/
      tid = omp_get_thread_num();
      if (tid == 0)
	{
	  nthreads = omp_get_num_threads();
	  Rprintf("Number of threads = %d\n", nthreads); // is this thread-safe?
	}
      Rprintf("Thread %d starting...\n",tid);
      
#pragma omp for schedule(dynamic,1)    
      for (i = 0; i < n; i++) {
	//Rprintf("id:%d, tid=%d\n", i, tid);
	Sim sim;
	Person person = Person(i);
	sim.create_process(&person);
	sim.run_simulation();
	sim.clear();
      }
    }
    
    // output
    return wrap(true);
  }

  RcppExport SEXP callSimTest2(SEXP parmsIn) {
    
    List parms(parmsIn);
    int n = as<int>(parms["n"]);
    
    Sim sim;
    Person person;

    for (int i = 0; i < n; i++) {
      person = Person(i);
      sim.create_process(&person);
      sim.run_simulation();
      sim.clear();
    }
    
    // output
    return wrap(true);
  }

  
} // anonymous namespace


/*
require(microsimulation)
.Call("callSimTest",list(n=10),package="microsimulation")

 */
