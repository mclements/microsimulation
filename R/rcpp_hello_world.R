#' microsimulation
#'
#' Discrete event simulations in both R and C++ with Tools for Cost-Effectiveness Analysis.
#'
#' @section Introduction:
#'
#' Discrete event simulations in both R and C++ with Tools for Cost-Effectiveness Analysis.
#'
#' @name microsimulation-package
#' @aliases microsimulation
#' @author Mark Clements \email{mark.clements@ki.se}
#' @references \url{https://github.com/mclements/microsimulation}
#' @seealso \code{\link[Rcpp]{sourceCpp}}
#' @useDynLib microsimulation, .registration=TRUE
#' @import Rcpp
#' @import methods
#' @importFrom graphics lines plot
#' @importFrom stats predict rnorm sd
#' @importFrom ascii ascii
"_PACKAGE"


#' Cat a string for the library archive for use in loading the package
#' 
#' @rdname Utilities
#' @return No return value, called for side effects
#' @export
LdFlags <- function()
    cat(system.file("lib/libmicrosimulation.a",  package="microsimulation", mustWork=TRUE))


#' On entry to the package, initialise the current stream in C.
#'
#' Is this function needed? We could define the current stream in open code.
#'
#' @param PACKAGE package from which this is called.
#'
#' @rdname Utilities
#' @return No return value, called for side effects
#' @export
microsimulation.init <- function(PACKAGE="microsimulation") {
    .C("r_create_current_stream",PACKAGE=PACKAGE)
    return(1)
}

#' On exit from the package, remove the current stream.
#'
#' Again, is this needed?
#' 
#' @param PACKAGE package from which this is called.
#'
#' @return No return value, called for side effects
#' @rdname Utilities
#' @export
microsimulation.exit <- function(PACKAGE="microsimulation") {
  .C("r_remove_current_stream",PACKAGE=PACKAGE)
  return(1)
}

## http://r.789695.n4.nabble.com/How-to-construct-a-valid-seed-for-l-Ecuyer-s-method-with-given-Random-seed-td4656340.html
#' Convert from signed to unsigned
#'
#' @param seed signed seed (possibly a vector)
#' @return unsigned seed
#' @rdname Utilities
#' @export 
unsigned <- function(seed) ifelse(seed < 0, seed + 2^32, seed)

#' Convert from unsigned to signed
#'
#' @param seed unsigned seed (possibly a vector)
#' @return signed seed
#' 
#' @rdname Utilities
#' @export 
signed <- function(seed) ifelse(seed>2^31, seed-2^32, seed)

#' Random draw for a positive (or other lower bound) random normal distribution
#'
#' @param n integer for the number of draws
#' @param mean numeric for the mean of the (untruncated) normal distribution (default=0)
#' @param sd numeric for the sd of the (untruncated) normal distribution (default=1)
#' @param lbound numeric for the lower bound (default=0)
#' @return numeric vector
#'
#' @importFrom stats rnorm sd
#' @rdname Utilities
#' @export
rnormPos <- function(n,mean=0,sd=1,lbound=0) {
    if (length(mean)<n) mean <- rep(mean,length=n)
    if (length(sd)<n) sd <- rep(sd,length=n)
    x <- rnorm(n,mean,sd)
    while(any(i <- which(x<lbound)))
        x[i] <- rnorm(length(i),mean[i],sd[i])
    x
}

#' Set the RngStream random number seed
#'
#' @param seed random number seed
#' @param PACKAGE package for the seed
#' @return invisibly returns the new seed
#' @rdname Utilities
#' @export
set.user.Random.seed <- function (seed,PACKAGE="microsimulation") {
  seed <- as.double(unsigned(seed))
  if (length(seed) == 1) seed <- rep(seed,6)
  if (length(seed) == 7) seed <- seed[-1]
  .C("r_set_user_random_seed",seed = seed,PACKAGE=PACKAGE)
  return(invisible(seed))
}

#' Advance the RngStream random number seed
#'
#' @param seed random number seed
#' @param n number of sub-streams to advance
#' @param PACKAGE package for the seed
#' @return the advanced seed
#'
#' @rdname Utilities
#' @export
advance.substream <- function (seed,n,PACKAGE="microsimulation") {
  seed <- as.double(unsigned(seed))
  if (length(seed) == 1) seed <- rep(seed,6)
  if (length(seed) == 7) seed <- seed[-1]
  .C("r_rng_advance_substream", seed = seed, n = as.integer(n), PACKAGE=PACKAGE)$seed
}

#' Advance the RngStream random number seed to the next sub-stream
#'
#' @param seed random number seed
#' @param PACKAGE package for the seed
#' @return invisibly returns TRUE -- called for side effect
#'
#' @rdname Utilities
#' @export
next.user.Random.substream <- function(PACKAGE="microsimulation") {
  .C("r_next_rng_substream", PACKAGE=PACKAGE)
  return(invisible(TRUE))
}

#' Get the current RngStream random seed
#'
#' @param PACKAGE package for the seed
#' @return random seed
#' @rdname Utilities
#' @export
user.Random.seed <- function(PACKAGE="microsimulation") {
  c(407L,
    as.integer(signed(.C("r_get_user_random_seed", seed=as.double(rep(1,6)),
                         PACKAGE=PACKAGE)$seed)))
}

#' Add labels to a (default) zero-based integer to give a factor
#'
#' @param obj integer or logical for factor levels
#' @param labels labels for the factor levels
#' @param start first value of the levels
#' @return the new factor
#' 
#' @rdname Utilities
#' @export
enum <- function(obj, labels, start=0) {
    ## stopifnot(is.logical(obj) || is.numeric(obj) || is.factor(obj))
    obj <- as.integer(obj)
    structure(factor(obj, levels=start + (0:(length(labels)-1)), labels=labels),start=start)
}

#' Assign labels to a (default) zero-based integer to give a factor
#'
#' @param obj integer or logical for factor levels
#' @param value labels for the factor levels
#' @return update the factor
#' @rdname Utilities
#' @export
"enum<-" <- function(obj, value) {
  enum(if(is.factor(obj)) unclass(obj)-1+attr(obj,"start") else obj, value)
}

#' S3 class for maintain state for .Random.seed
#'
#' @return a list with oldseed (the old value of .Random.seed), and reset(), which resets .Random.seed
#'
#' @rdname Utilities
#' @export
RNGstate <- function() {
  ## house-keeping for random streams (see parallel::clusterSetRNGStream)
  oldseed <- if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else NULL
  oldkind <- RNGkind()
  reset <- function() {
    if (!is.null(oldseed))
      assign(".Random.seed", oldseed, envir = .GlobalEnv)
    else {
        ## clean up if we created a .Random.seed
        if (exists(".Random.seed" ,envir = .GlobalEnv, inherits = FALSE))
            rm(.Random.seed, envir = .GlobalEnv, inherits = FALSE)
            do.call(RNGkind, as.list(oldkind))
    }
  }
  list(oldseed = oldseed, reset = reset)
}

#' Utility to calculate the cost-efficiency frontier
#'
#' @param x vector of x coordinates
#' @param y vector y coordinates
#' @param concave logical for whether to calculate a concave frontier (default=TRUE)
#' @param convex logical for whether to calculate a convex frontier (default=NULL)
#' @return a list with components x and y for the frontier
#' @rdname Utilities
#' @export
frontier <- function (x, y, concave=TRUE, convex=NULL)
{
    ## check arguments
    stopifnot(is.logical(concave),
              is.null(convex) || is.logical(convex),
              is.numeric(x),
              is.numeric(y),
              length(x) == length(y),
              !any(is.na(x)),
              !any(is.na(y)))
    ## Change concave if convex is defined
    if (!is.null(convex))
        concave <- !convex
    if (concave) {
        ichull <- grDevices::chull(cbind(x, y))
        ## case: length == 0
        if ((n <- length(ichull)) == 0)
            return(NULL)
        ## case: length == 1
        if (n == 1)
            return(list(x=x[ichull], y=y[ichull]))
        ## case: length > 1
        ## if min(x) value is not first in ichull, then re-order
        if ((iminx <- which.min(x[ichull])) > 1) {
            ichull <- ichull[c(iminx:n, 1:(iminx-1))]
        }
        ## find those that are increasing for x and y
        include <- c(TRUE,diff(x[ichull])>0 & diff(y[ichull])>0)
        list(x=x[ichull][include], y=y[ichull][include])
    } else # convex case
        with(Recall(y,x,concave=TRUE), list(x=y, y=x))
}

#' plot lines for a frontier
#'
#' @param x vector of x coordinates
#' @param y vector of y coordinates
#' @param pch type of pch for the plotted symbols (default=19)
#' @param type join type (default="b")
#' @param ... other arguments to lines
#' 
#' @return No return value, called for side effects
#' 
#' @rdname Utilities
lines_frontier <- function(x,y,pch=19,type="b",...) {
    index <- frontier(x,y)
    lines(x[index],y[index],pch=pch,type=type,...)
}

#' Integrate a discounted value
#'
#' @param y the undiscounted value
#' @param start the start time
#' @param finish the finish time
#' @param dr discount rate, expressed as a percentage
#' @return numeric discounted value
#' @export
discountedInterval <- function(y, start, finish, dr) {
    duration = finish - start
    ifelse(dr<=0 | duration==0, y*duration,
           y / (1 + dr)^ start / log(1 + dr) *
           (1 - 1 / (1 + dr)^ duration))
}

#' Discounted value
#'
#' @param y the undiscounted value
#' @param time the time of the event
#' @param dr discount rate, expressed as a percentage
#' @return numeric vector
#' @rdname Utilities
#' @export
discountedPoint <- function(y, time, dr)
    ifelse(dr<=0 | time==0, y, y / (1 + dr) ^ time)


#' Call the Person simulation from C++
#'
#' Example that uses the RngStream random number generator
#' 
#' @param n number of simulations (default=20)
#' @param seed random number seed
#' @return data-frame
#' @rdname Examples
#' @export
callPersonSimulation <- function(n=20,seed=rep(12345,6)) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("user")
  set.user.Random.seed(seed)
  stateT =c("Healthy","Localised","DxLocalised","LocallyAdvanced",
    "DxLocallyAdvanced","Metastatic","DxMetastatic","Death")
  eventT =c("toDeath", "toPCDeath", "toLocalised", "toDxLocalised",
	      "toDxLocallyAdvanced",
	      "toLocallyAdvanced", "toMetastatic", "toDxMetastatic")
  out <- .Call(.callPersonSimulation,
               as.integer(rep(seed,length=6)), # magic
               list(n=as.integer(n)))
  out$state = enum(out$state,stateT)
  out$event = enum(out$event,eventT)
  as.data.frame(out)
}

#' Call the SimplePerson simulation from C++
#'
#' Example that uses the Mersenne-Twister random number generator
#' 
#' @param n number of simulations (default=10)
#' @return data-frame
#' @rdname Examples
#' @export
callSimplePerson <- function(n=10) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  set.seed(12345)
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call(.callSimplePerson,
               parms=list(n=as.integer(n)))
  out$state = enum(out$state,stateT)
  out$event = enum(out$event,eventT)
  as.data.frame(out)
}

#' Call the SimplePerson2 simulation from C++
#'
#' Example that uses the Mersenne-Twister random number generator
#' 
#' @param n number of simulations (default=10)
#' @return data-frame
#' @rdname Examples
#' @export
callSimplePerson2 <- function(n=10) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  set.seed(12345)
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call(.callSimplePerson2,
               parms=list(n=as.integer(n)))
  reader <- function(obj)
      cbind(data.frame(state=enum(obj[[1]],stateT)),
          data.frame(obj[-1]))
  out <- lapply(out,reader)
  out$events <- with(out$events, data.frame(state=state,event=enum(event,eventT),age=age,number=number))
  out$pt <- with(out$pt, data.frame(state=state,age=age,pt=pt))
  out$prev <- with(out$prev, data.frame(state=state,age=age,number=number))
  out
}

#' Call the IllnessDeath simulation from C++
#'
#' Example that uses the Mersenne-Twister random number generator
#' 
#' @param n number of simulations (default=10)
#' @param cure probability of cure
#' @param zsd frailty standard deviation
#' @return data-frame
#' @rdname Examples
#' @export
callIllnessDeath <- function(n=10L,cure=0.1,zsd=0) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  set.seed(12345)
  stateT <- c("Healthy","Cancer")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call(.callIllnessDeath,
               parms=list(n=as.integer(n),cure=as.double(cure),zsd=as.double(zsd)))
  reader <- function(obj)
      cbind(data.frame(state=enum(obj[[1]],stateT)),
          data.frame(obj[-1]))
  out <- lapply(out,reader)
  out$events <- with(out$events, data.frame(state=state,event=enum(event,eventT),age=age,number=number))
  out$pt <- with(out$pt, data.frame(state=state,age=age,pt=pt))
  out$prev <- with(out$prev, data.frame(state=state,age=age,number=number))
  out
}

#' Calculate incremental cost-effectiveness ratios from two objects
#'
#' @param object1 first object
#' @param object2 second object
#' @param ... other arguments
#' @rdname Utilities
#' @export 
ICER <- function(object1, object2, ...)
    UseMethod("ICER")

#' Reference class implementation of an event queue in R
#'
#'   This event queue is simple and useful for pedagogic purposes.
#' 
#'   The algorithm for pushing values into the queue is computationally
#'   very simple: simply rank the times using \code{order()} and re-order
#'   times and events. This approach is probably of acceptable performance
#'   for smaller queue. A more computationally efficient approach for
#'   pushing into larger queues would be to use a binary search (e.g. using
#'   \code{findInterval()}).
#' 
#' For faster alternatives, see \code{pqueue} and \code{PQueueRef}.
#'
#' @examples
#' pq = new("EventQueue")
#' pq$push(3,"Clear drains")
#' pq$push(4, "Feed cat")
#' pq$push(5, "Make tea")
#' pq$push(1, "Solve RC tasks")
#' pq$push(2, "Tax return")
#' while(!pq$empty())
#'   print(pq$pop())
#'
#' @import methods
#' @exportClass EventQueue
#' @field times vector of times
#' @field events list of events
#' @rdname Classes
EventQueue <-
    setRefClass("EventQueue",
                fields = list(times = "numeric", events = "list"),
                methods = list(
                    help = function() {
                        'Reference class implementation of an event queue. Fields for the event times and the events list.'
                    },
                    push = function(time,event) {
                        'Method to insert the event at the given time'
                        insert.ord <- findInterval(time,times)
                        times <<- append(times,time,insert.ord)
                        events <<- append(events,list(event),insert.ord)
                    },
                    pop = function() {
                        'Method to remove the head of the event queue and return its value'
                        head <- list(event=events[[1]], time=times[1])
                        times <<- times[-1]
                        events <<- events[-1]
                        return(head)
                    },
                    empty = function() {
                        'Method to check whether there are no events in the queue'
                        length(times) == 0
                        },
                    clear = function() {
                        'Method to clear the event queue'
                        times <<- numeric()
                        events <<- list()
                    },
                    cancel = function(predicate, ...) {
                        'Method to remove events that satisfy some predicate'
                        if (!empty()) {
                            i <- sapply(events, predicate, ...)
                            stopifnot(is.logical(i))
                            i[is.na(i)] <- TRUE
                            times <<- times[!i]
                            events <<- events[!i]
                        }
                    }))

#' Reference class implementation of a discrete event simulation
#'
#' Inherit from this class to represent a discrete event simulation. The
#' API is similar to that for Omnet++, where an \code{init} method sets up
#' the initial events using the \code{scheduleAt(time,event)} method, the
#' messages are handled using the \code{handleMessage(event)} method, the
#' simulation is run using the \code{run} method, and the \code{final}
#' method is called at the end of the simulation.
#' 
#' @examples
#' DES = setRefClass("DES",
#'                   contains = "BaseDiscreteEventSimulation",
#'                   methods=list(
#'                       init=function() {
#'                          scheduleAt(3,"Clear drains")
#'                          scheduleAt(4, "Feed cat")
#'                          scheduleAt(5, "Make tea")
#'                          scheduleAt(1, "Solve RC tasks")
#'                          scheduleAt(2, "Tax return")
#'                       },
#'                       handleMessage=function(event) print(event)))
#'
#' des = new("DES")
#' des$run()
#' @examples
#' \dontrun{
#' testRsimulation1 <- function() {
#'     ## A simple example
#'     Simulation <-
#'         setRefClass("Simulation",
#'                     contains = "BaseDiscreteEventSimulation")
#'     Simulation$methods(
#'         init = function() {
#'             scheduleAt(rweibull(1,8,85), "Death due to other causes")
#'             scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
#'         },
#'         handleMessage = function(event) {
#'             if (event %in% c("Death due to other causes", "Cancer death")) {
#'                 clear()
#'                 print(event)
#'             }
#'             else if (event == "Cancer diagnosis") {
#'                 if (runif(1) < 0.5)
#'                     scheduleAt(now() + rweibull(1,2,10), "Cancer death")
#'                 print(event)
#'             }
#'         })
#'     Simulation$new()$run()
#' }
#' 
#' ## An extension with individual life histories
#' testRsimulation2 <- function(n=100) {
#'     Simulation <-
#'         setRefClass("Simulation",
#'                     contains = "BaseDiscreteEventSimulation",
#'                     fields = list(state = "character", report = "data.frame"))
#'     Simulation$methods(
#'         init = function() {
#'             report <<- data.frame()
#'             state <<- "Healthy"
#'             scheduleAt(rweibull(1,8,85), "Death due to other causes")
#'             scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
#'         },
#'         handleMessage = function(event) {
#'             report <<- rbind(report, data.frame(state = state,
#'                                                 begin = attr(event,"sendingTime"),
#'                                                 end = currentTime,
#'                                                 event = event,
#'                                                 stringsAsFactors = FALSE))
#'             if (event %in% c("Death due to other causes", "Cancer death")) {
#'                 clear()
#'             }
#'             else if (event == "Cancer diagnosis") {
#'                 state <<- "Cancer"
#'                 if (runif(1) < 0.5)
#'                     scheduleAt(now() + rweibull(1,2,10), "Cancer death")
#'             }
#'         },
#'         final = function() report)
#'     sim <- Simulation$new()
#'     do.call("rbind", lapply(1:n, function(id) data.frame(id=id,sim$run())))
#' }
#' 
#' ## reversible illness-death model
#' testRsimulation3 <- function(n=100) {
#'     Simulation <-
#'         setRefClass("Simulation",
#'                     contains = "BaseDiscreteEventSimulation",
#'                     fields = list(state = "character", everCancer = "logical",
#'                                   report = "data.frame"))
#'     Simulation$methods(
#'         init = function() {
#'             report <<- data.frame()
#'             state <<- "Healthy"
#'             everCancer <<- FALSE
#'             scheduleAt(rweibull(1,8,85), "Death due to other causes")
#'             scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
#'         },
#'         handleMessage = function(event) {
#'             report <<- rbind(report, data.frame(state = state,
#'                                                 everCancer = everCancer,
#'                                                 begin = attr(event,"sendingTime"),
#'                                                 end = currentTime,
#'                                                 event = event,
#'                                                 stringsAsFactors = FALSE))
#'             if (event %in% c("Death due to other causes", "Cancer death")) {
#'                 clear()
#'             }
#'             else if (event == "Cancer diagnosis") {
#'                 state <<- "Cancer"
#'                 everCancer <<- TRUE
#'                 if (runif(1) < 0.5)
#'                     scheduleAt(now() + rweibull(1,2,10), "Cancer death")
#'                 scheduleAt(now() + 10, "Recovery")
#'             }
#'             else if (event == "Recovery") {
#'                 state <<- "Healthy"
#'                 scheduleAt(now() + rexp(1,10), "Cancer diagnosis")
#'             }
#'         },
#'         final = function() report)
#'     sim <- Simulation$new()
#'     do.call("rbind", lapply(1:n, function(id) data.frame(id=id,sim$run())))
#' }
#' 
#' ## cancer screening
#' testRsimulation4 <- function(n=1) {
#'     Simulation <-
#'         setRefClass("Simulation",
#'                     contains = "BaseDiscreteEventSimulation",
#'                     fields = list(state = "character", report = "data.frame"))
#'     Simulation$methods(
#'         init = function() {
#'             report <<- data.frame()
#'             state <<- "Healthy"
#'             scheduleAt(rweibull(1,8,85), "Death due to other causes")
#'             scheduleAt(rweibull(1,3,90), "Cancer onset")
#'             scheduleAt(50,"Screening")
#'         },
#'         handleMessage = function(event) {
#'             report <<- rbind(report, data.frame(state = state,
#'                                                 begin = attr(event,"sendingTime"),
#'                                                 end = currentTime,
#'                                                 event = event,
#'                                                 stringsAsFactors = FALSE))
#'             if (event %in% c("Death due to other causes", "Cancer death")) {
#'                 clear()
#'             }
#'             else if (event == "Cancer onset") {
#'                 state <<- event
#'                 dx <- now() + rweibull(1,2,10)
#'                 scheduleAt(dx, "Clinical cancer diagnosis")
#'                 scheduleAt(dx + rweibull(1,1,10), "Cancer death")
#'                 scheduleAt(now() + rweibull(1,1,10), "Metastatic cancer")
#'             }
#'             else if (event == "Metastatic cancer") {
#'                 state <<- event
#'                 cancel(function(event) event %in%
#'                        c("Clinical cancer diagnosis","Cancer death")) # competing events
#'                 scheduleAt(now() + rweibull(1,2,5), "Cancer death")
#'             }
#'             else if (event == "Clinical cancer diagnosis") {
#'                 state <<- event
#'                 cancel(function(event) event == "Metastatic cancer")
#'             }
#'             else if (event == "Screening") {
#'                 switch(state,
#'                        "Cancer onset" = {
#'                            state <<- "Screen-detected cancer diagnosis"
#'                            cancel(function(event) event %in%
#'                                   c("Clinical cancer diagnosis","Metastatic cancer"))
#'                        },
#'                        "Metastatic cancer" = {}, # ignore
#'                        "Clincal cancer diagnosis" = {}, # ignore
#'                        "Healthy" = {
#'                            if (now()<=68) scheduleAt(now()+2, "Screening")
#'                        })
#'             }
#'             else stop(event)
#'         },
#'         final = function() report)
#'     sim <- Simulation$new()
#'     do.call("rbind", lapply(1:n, function(id) data.frame(id=id,sim$run())))
#' }
#' 
#' ## ticking bomb - toy example
#' testRsimulation5 <- function(n=1) {
#'     Simulation <-
#'         setRefClass("Simulation",
#'                     contains = "BaseDiscreteEventSimulation",
#'                     fields = list(report = "data.frame"))
#'     Simulation$methods(
#'         init = function() {
#'             report <<- data.frame()
#'             scheduleAt(rexp(1,1), "tick")
#'             if (runif(1)<0.1)
#'                 scheduleAt(rexp(1,1), "explosion")
#'         },
#'         handleMessage = function(event) {
#'             report <<- rbind(report, data.frame(begin = attr(event,"sendingTime"),
#'                                                 end = currentTime,
#'                                                 event = event,
#'                                                 stringsAsFactors = FALSE))
#'             if (event == "explosion")
#'                 clear()
#'             else {
#'                 clear() # queue
#'                 if (event == "tick") scheduleAt(currentTime+rexp(1,1), "tock")
#'                 else scheduleAt(currentTime+rexp(1,1), "tick")
#'                 if (runif(1)<0.1)
#'                     scheduleAt(currentTime+rexp(1,1), "explosion")
#'             }
#'         },
#'         final = function() report)
#'     sim <- Simulation$new()
#'     do.call("rbind", lapply(1:n, function(id) data.frame(id=id,sim$run())))
#' }
#' }
#' 
#' @import methods
#' @exportClass EventQueue
#' @field times vector of times
#' @field events list of events
#' @rdname Classes
BaseDiscreteEventSimulation <-
    setRefClass("BaseDiscreteEventSimulation",
                contains = "EventQueue",
                ## contains = "PQueueRef",
                fields = list(currentTime = "numeric",
                    previousEventTime = "numeric"),
                methods = list(
                    help = function() {
                        'Reference class implementation of an event-oriented discrete event simulation. Fields for the event times and the events list.'
                    },
                    scheduleAt = function(time, event) {
                        'Method that adds attributes for the event time and the sendingTime, and then insert the event into the event queue'
                        attr(event,"time") <- time
                        attr(event,"sendingTime") <- currentTime
                        push(time, event)
                    },
                    init = function() {
                        'Virtual method to initialise the event queue and attributes'
                        stop("VIRTUAL!")
                    },
                    handleMessage = function(event) {
                        'Virtual method to handle the messages as they arrive'
                        stop("VIRTUAL!")
                    },
                    final = function() {
                        'Method for finalising the simulation'
                        NULL
                    },
                    now = function() currentTime,
                    previous = function() previousEventTime,
                    reset = function(startTime = 0.0) {
                        'Method to reset the event queue'
                        clear()
                        previousEventTime <<- currentTime <<- startTime
                    },
                    run = function(startTime = 0.0) {
                        'Method to run the simulation'
                        reset(startTime)
                        init()
                        while(!empty()) {
                            event <- pop()$event
                            currentTime <<- attr(event, "time")
                            handleMessage(event)
                            previousEventTime <<- currentTime
                        }
                        final()
                    }))

#' S3 class to work with RngStream objects
#'
#' @param nextStream whether to move to the next stream (default=TRUE)
#' @param iseed set seed after changing RNG (otherwise keep the current seed)
#' @return list of class \code{RNGStream} with components:
#' \describe{
#' \item{resetRNGkind}{function to reset to the previous RNG and seed}
#' \item{seed}{function to return the current seed}
#' \item{open}{function to use the current seed}
#' \item{close}{function to make the current seed equal to .Random.seed}
#' \item{resetStream}{function to move back to start of stream}
#' \item{resetSubStream}{function to move back to start of sub-stream}
#' \item{nextSubStream}{function to move to next sub-stream}
#' \item{nextStream}{function to move to next stream}
#' }
#' @examples
#' ## set up one stream
#' s1 <- RNGStream()
#' s1$open()
#' rnorm(1)
#' s1$nextSubStream()
#' rnorm(1)
#' ## reset the stream
#' s1$resetStream()
#' rnorm(2)
#' s1$nextSubStream()
#' rnorm(2)
#' 
#' ## now do with two streams
#' s1$resetStream()
#' s2 <- RNGStream()
#' with(s1,rnorm(1))
#' with(s2,rnorm(1))
#' s1$nextSubStream()
#' with(s1,rnorm(1))
#' ## now reset the streams and take two samples each time
#' s1$resetStream()
#' s2$resetStream()
#' with(s1,rnorm(2))
#' with(s2,rnorm(2))
#' s1$nextSubStream()
#' with(s1,rnorm(2))
#' @rdname RNGStream
#' @export
RNGStream <- function(nextStream = TRUE, iseed = NULL) {
  stopifnot(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) || !is.null(iseed))
  oldseed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else NULL
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  if (!is.null(iseed)) set.seed(iseed)
  current <- if (nextStream) parallel::nextRNGStream(.Random.seed) else .Random.seed
  .Random.seed <<- startOfStream <- startOfSubStream <- current
  structure(list(resetRNGkind = function()
  {
      if (!is.null(oldseed))
          assign(".Random.seed", oldseed, envir = .GlobalEnv)
      else   rm(.Random.seed, envir = .GlobalEnv)
  },
  seed = function() current,
  open = function() .Random.seed <<- current,
  close = function() current <<- .Random.seed,
  resetStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream,
  resetSubStream = function() .Random.seed <<- current <<- startOfSubStream,
  nextSubStream = function() .Random.seed <<- current <<- startOfSubStream <<- parallel::nextRNGSubStream(startOfSubStream),
  nextStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream <<- parallel::nextRNGStream(startOfStream)),
  class="RNGStream")
}

#' Use RNGStream as an old class
#'
#' @name RNGStream-class
#' @aliases RNGStream
#' @rdname RNGStream
#' @exportClass RNGStream
setOldClass("RNGStream") 

#' With method for RNGStream S3 class
#'
#' @param data object of type RNGStream
#' @param expr expression using the RNGStream
#' @param ... other arguments passed to eval()
#' @return the value from the expression
#' @rdname RNGStream
with.RNGStream <- function(data,expr,...) {
  data$open()
  out <- eval(substitute(expr), enclos = parent.frame(), ...)
  data$close()
  out
}

#' Old data used in the prostata model
#'@rdname Data
"fhcrcData"

#' C++ function
#' @rdname Internal
#' @name callCalibrationSimulation
#' @return data-frame
NULL

#' C++ function
#' @rdname Internal
#' @name r_create_current_stream
#' @return No return value, called for side effects
#' @export
NULL

#' C++ function
#' @rdname Internal
#' @name r_remove_current_stream
#' @return No return value, called for side effects
#' @export
NULL

#' C++ function
#' @rdname Internal
#' @name r_set_user_random_seed
#' @return No return value, called for side effects
#' @export
NULL

#' C++ function
#' @rdname Internal
#' @name r_rng_advance_substream
#' @return No return value, called for side effects
#' @export
NULL

#' C++ function
#' @rdname Internal
#' @name r_next_rng_substream
#' @return No return value, called for side effects
#' @export
NULL

#' C++ function
#' @rdname Internal
#' @name r_get_user_random_seed
#' @return No return value, called for side effects
#' @export
NULL

