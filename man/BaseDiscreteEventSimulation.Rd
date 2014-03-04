\name{BaseDiscreteEventSimulation}
\docType{methods}
\alias{class-BaseDiscreteEventSimulation}
\alias{new,BaseDiscreteEventSimulation-method}
\alias{scheduleAt,BaseDiscreteEventSimulation-method}
\alias{init,BaseDiscreteEventSimulation-method}
\alias{handleMessage,BaseDiscreteEventSimulation-method}
\alias{final,BaseDiscreteEventSimulation-method}
\alias{now,BaseDiscreteEventSimulation-method}
\alias{reset,BaseDiscreteEventSimulation-method}
\alias{run,BaseDiscreteEventSimulation-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An R5 reference class to implement discrete event simulation. 
}
\description{
Inherit from this class to represent a discrete event simulation. The
API is similar to that for Omnet++, where an \code{init} method sets up
the initial events using the \code{scheduleAt(time,event)} method, the
messages are handled using the \code{handleMessage(event)} method, the
simulation is run using the \code{run} method, and the \code{final}
method is called at the end of the simulation.
}
\usage{
\S4method{new}{BaseDiscreteEventSimulation}(currentTime,previousEventTime)
\S4method{scheduleAt}{BaseDiscreteEventSimulation}(time,event)
\S4method{init}{BaseDiscreteEventSimulation}()
\S4method{handleMessage}{BaseDiscreteEventSimulation}(event)
\S4method{final}{BaseDiscreteEventSimulation}()
\S4method{now}{BaseDiscreteEventSimulation}()
\S4method{reset}{BaseDiscreteEventSimulation}(startTime = 0.0)
\S4method{run}{BaseDiscreteEventSimulation}(startTime = 0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{currentTime}{
    Simulation-level variable recording the current time. This is also
    returned by the \code{now} method.
  }
  \item{previousEventTime}{Simulation-level variable recording the
    previous event time.}
  \item{time}{schedule time that the event will take place}
  \item{event}{A value that represents the type of event. How this value
    is used is determined by the \code{handleMessage} method.}
  \item{startTime}{From when the simulation starts (or is reset)}
}
\details{
A concrete discrete event simulation class would inherit from this
class; see below for examples. The methods for \code{init} and
\code{handleMessage} are required to be implemented in the concrete class.

The \code{scheduleAt} method adds attributes to the event for the
\code{time} for the scheduled event and the \code{sendingTime} when the
event was originally scheduled. 

The implementation inherits from the R5 reference class
\code{EventQueue}. The events are inserted using an ordered vector
(rather than something more sophisticated, such as \code{findInterval}).

}
%% \value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%% }
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~
%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
\dontrun{
testRsimulation1 <- function() {
    ## A simple example
    Simulation <-
        setRefClass("Simulation",
                    contains = "BaseDiscreteEventSimulation")
    Simulation$methods(
        init = function() {
            scheduleAt(rweibull(1,8,85), "Death due to other causes")
            scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
        },
        handleMessage = function(event) {
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
    Simulation$new()$run()
}

## An extension with individual life histories
testRsimulation2 <- function(n=100) {
    Simulation <-
        setRefClass("Simulation",
                    contains = "BaseDiscreteEventSimulation",
                    fields = list(state = "character", report = "data.frame"))
    Simulation$methods(
        init = function() {
            report <<- data.frame()
            state <<- "Healthy"
            scheduleAt(rweibull(1,8,85), "Death due to other causes")
            scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
        },
        handleMessage = function(event) {
            report <<- rbind(report, data.frame(state = state,
                                                begin = attr(event,"sendingTime"),
                                                end = currentTime,
                                                event = event,
                                                stringsAsFactors = FALSE))
            if (event \%in\% c("Death due to other causes", "Cancer death")) {
                clear()
            }
            else if (event == "Cancer diagnosis") {
                state <<- "Cancer"
                if (runif(1) < 0.5)
                    scheduleAt(now() + rweibull(1,2,10), "Cancer death")
            }
        },
        final = function() report)
    sim <- Simulation$new()
    do.call("rbind", lapply(1:n, function(id) data.frame(id=id,sim$run())))
}

## reversible illness-death model
testRsimulation3 <- function(n=100) {
    Simulation <-
        setRefClass("Simulation",
                    contains = "BaseDiscreteEventSimulation",
                    fields = list(state = "character", everCancer = "logical", report = "data.frame"))
    Simulation$methods(
        init = function() {
            report <<- data.frame()
            state <<- "Healthy"
            everCancer <<- FALSE
            scheduleAt(rweibull(1,8,85), "Death due to other causes")
            scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
        },
        handleMessage = function(event) {
            report <<- rbind(report, data.frame(state = state,
                                                everCancer = everCancer,
                                                begin = attr(event,"sendingTime"),
                                                end = currentTime,
                                                event = event,
                                                stringsAsFactors = FALSE))
            if (event \%in\% c("Death due to other causes", "Cancer death")) {
                clear()
            }
            else if (event == "Cancer diagnosis") {
                state <<- "Cancer"
                everCancer <<- TRUE
                if (runif(1) < 0.5)
                    scheduleAt(now() + rweibull(1,2,10), "Cancer death")
                scheduleAt(now() + 10, "Recovery")
            }
            else if (event == "Recovery") {
                state <<- "Healthy"
                scheduleAt(now() + rexp(1,10), "Cancer diagnosis")
            }
        },
        final = function() report)
    sim <- Simulation$new()
    do.call("rbind", lapply(1:n, function(id) data.frame(id=id,sim$run())))
}

## cancer screening
testRsimulation4 <- function(n=1) {
    Simulation <-
        setRefClass("Simulation",
                    contains = "BaseDiscreteEventSimulation",
                    fields = list(state = "character", report = "data.frame"))
    Simulation$methods(
        init = function() {
            report <<- data.frame()
            state <<- "Healthy"
            scheduleAt(rweibull(1,8,85), "Death due to other causes")
            scheduleAt(rweibull(1,3,90), "Cancer onset")
            scheduleAt(50,"Screening")
        },
        handleMessage = function(event) {
            report <<- rbind(report, data.frame(state = state,
                                                begin = attr(event,"sendingTime"),
                                                end = currentTime,
                                                event = event,
                                                stringsAsFactors = FALSE))
            if (event \%in\% c("Death due to other causes", "Cancer death")) {
                clear()
            }
            else if (event == "Cancer onset") {
                state <<- event
                dx <- now() + rweibull(1,2,10)
                scheduleAt(dx, "Clinical cancer diagnosis")
                scheduleAt(dx + rweibull(1,1,10), "Cancer death")
                scheduleAt(now() + rweibull(1,1,10), "Metastatic cancer")
            }
            else if (event == "Metastatic cancer") {
                state <<- event
                remove(function(event) event %in% c("Clinical cancer diagnosis","Cancer death")) # competing events
                scheduleAt(now() + rweibull(1,2,5), "Cancer death")
            }
            else if (event == "Clinical cancer diagnosis") {
                state <<- event
                remove(function(event) event == "Metastatic cancer")
            }
            else if (event == "Screening") {
                switch(state,
                       "Cancer onset" = {
                           state <<- "Screen-detected cancer diagnosis"
                           remove(function(event) event %in% c("Clinical cancer diagnosis","Metastatic cancer"))
                       },
                       "Metastatic cancer" = {}, # ignore
                       "Clincal cancer diagnosis" = {}, # ignore
                       "Healthy" = {
                           if (now()<=68) scheduleAt(now()+2, "Screening")
                       })
            }
            else stop(event)
        },
        final = function() report)
    sim <- Simulation$new()
    do.call("rbind", lapply(1:n, function(id) data.frame(id=id,sim$run())))
}

## ticking bomb - toy example
testRsimulation5 <- function(n=1) {
    Simulation <-
        setRefClass("Simulation",
                    contains = "BaseDiscreteEventSimulation",
                    fields = list(report = "data.frame"))
    Simulation$methods(
        init = function() {
            report <<- data.frame()
            scheduleAt(rexp(1,1), "tick")
            if (runif(1)<0.1)
                scheduleAt(rexp(1,1), "explosion")
        },
        handleMessage = function(event) {
            report <<- rbind(report, data.frame(begin = attr(event,"sendingTime"),
                                                end = currentTime,
                                                event = event,
                                                stringsAsFactors = FALSE))
            if (event == "explosion")
                clear()
            else {
                clear() # queue
                if (event == "tick") scheduleAt(currentTime+rexp(1,1), "tock")
                else scheduleAt(currentTime+rexp(1,1), "tick")
                if (runif(1)<0.1)
                    scheduleAt(currentTime+rexp(1,1), "explosion")
            }
        },
        final = function() report)
    sim <- Simulation$new()
    do.call("rbind", lapply(1:n, function(id) data.frame(id=id,sim$run())))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
