\name{BaseDiscreteEventSimulation}
\docType{methods}
\alias{class-BaseDiscreteEventSimulation}
\alias{scheduleAt,BaseDiscreteEventSimulation-method}
\alias{init,BaseDiscreteEventSimulation-method}
\alias{handleMessage,BaseDiscreteEventSimulation-method}
\alias{now,BaseDiscreteEventSimulation-method}
\alias{run,BaseDiscreteEventSimulation-method}
\title{
  Base R5 class for discrete event simulation
}
\description{
  Base R5 class for discrete event simulation. It requires that an
  inherited R5 class defines the init()
  and handleMessage(event) methods (and possibly the final()
  method). See below for a simple example. 
}
\usage{
\S4method{scheduleAt}{BaseDiscreteEventSimulation}(time, event)
\S4method{init}{BaseDiscreteEventSimulation}()
\S4method{handleMessage}{BaseDiscreteEventSimulation}(event)
\S4method{now}{BaseDiscreteEventSimulation}()
\S4method{run}{BaseDiscreteEventSimulation}()
}
\arguments{
  \item{time}{Scheduled time for the event (numeric).}
  \item{event}{An arbitrary object representing an event (or
    message). The event will be handled by the \code{handleMessage(event)} method.
  }
}
\section{Methods}{
  \describe{
    \item{\code{init()}}{
      Initialise the simulation run, inserting one or more events using
      \code{scheduleAt(time,event).}
    }
    \item{\code{scheduleAt(time,event)}}{
      Insert an event at a given time into the event queue.
    }
    \item{\code{now()}}{
      Utility function to get the current simulation time.
    }
    \item{\code{handleMessage(event)}}{
      
    }
    \item{\code{run()}}{
    }
  }
  Other methods are inherited from \code{\link{class-EventQueue}}.
}
\details{
  
  This is a simple discrete event simulation system based on the
  event-oriented paradigm with a single process. An inherited class
  defines the \code{init()} and \code{handleMessage()} methods and then an instance can
  \code{run()}'s the model once or more. For microsimulation, the model would
  be run many times.
  
  The \code{now()} method provides the current time. There are two useful
  fields: \code{currentTime}, which is returned by \code{run()}, and
  \code{previousEventTime}.  The simulation model uses the \code{EventQueue}
  class. The \code{run()} method sets \code{previousEventTime}, calls \code{init()},
  then processes the event queue until it is empty, and finishes by
  calling \code{final()}. The inherited
  \code{clear()} and \code{remove(eventPredicate)} methods from the \code{EventQueue} class
  are also useful: the former clears the event queue, while the latter
  removes events from the queue that satisfy the event predicate.
  
  This is a lightweight class, where the form of the events (or
  "self-messages", in OMNeT++ parlance) is arbitrary and determined by
  the handleMessage() method. At present, there is little infrastructure
  for reporting from a simulation; any such report could be stored as a
  field in the inherited class and possibly returned by the \code{final()} method.
  
  For context, the event orientation is different to the process
  orientation found in \href{http://simpy.sourceforge.net/}{Simpy}. The API used here borrows from the \href{www.omnetpp.org}{OMNeT++}
  API for the
  \href{http://www.omnetpp.org/doc/omnetpp/api/classcSimpleModule.html}{cSimpleModule()}
  class (although we changed the OMNET++ \code{initialise()} and \code{finalise()}
  methods to \code{init()} and \code{final()} to avoid confusion with R5's
  \code{initialize()} method).
  
  This class is useful for introducing discrete event simulation and for
  prototyping, however it is not ideal for large-scale microsimulation
  due to slow execution time -- which largely explains the parallel
  development in C++.
}
% \value{
%   The output class can be
%   }
\examples{
\dontrun{
Simulation <-
  setRefClass("Simulation",
              contains = "BaseDiscreteEventSimulation")
Simulation$methods(init = function() {
  clear()
  scheduleAt(rweibull(1,8,85), "Death due to other causes")
  scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
})
Simulation$methods(handleMessage = function(event) {
  if (event \%in\% c("Death due to other causes", "Cancer death")) {
    clear()
    print(event)
  }
  else if (event == "Cancer diagnosis") {
    if (runif(1) < 0.5)
      scheduleAt(now() + rweibull(1,2,10), "Cancer death")
    print(event)
  }
})
set.seed(123)
sim <- Simulation$new()
sim$run()
sim$run()
}
}
\seealso{class-EventQueue}