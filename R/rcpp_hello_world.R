
rcpp_hello_world <- function(){
	.Call( "rcpp_hello_world", PACKAGE = "microsimulation" )
}

.microsimulation.init <- function () {
  .Call("r_create_current_stream",PACKAGE="microsimulation")
  return(1)
}

.microsimulation.exit <- function () {
  .Call("r_remove_current_stream",PACKAGE="microsimulation")
  return(1)
}

## callPersonSimulation <- function(n=500L)
##   .C("callPersonSimulation",as.integer(rep(12345,6)),
##      as.double(1),as.integer(n),
##      out=as.double(1:2),as.integer(2),PACKAGE="microsimulation")$out

enum <- function(obj, labels)
  factor(obj, levels=0:(length(labels)-1), labels=labels)

RNGstate <- function() {
  ## house-keeping for RNGStreams (see parallel::clusterSetRNGStream)
  oldseed <- if (exists(".Random.seed", envir = .GlobalEnv, 
                        inherits = FALSE)) 
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else NULL
  reset <- function() {
    if (!is.null(oldseed)) 
      assign(".Random.seed", oldseed, envir = .GlobalEnv)
    else rm(.Random.seed, envir = .GlobalEnv)
  }
  list(reset = reset)
}
  
callPersonSimulation <- function(n=20,seed=rep(12345,6)) {
  state <- RNGstate()
  RNGkind("user")
  stateT =c("Healthy","Localised","DxLocalised","LocallyAdvanced",
    "DxLocallyAdvanced","Metastatic","DxMetastatic","Death")
  eventT =c("toDeath", "toPCDeath", "toLocalised", "toDxLocalised",
	      "toDxLocallyAdvanced",
	      "toLocallyAdvanced", "toMetastatic", "toDxMetastatic")
  out <- .Call("callPersonSimulation",
               as.integer(rep(seed,length=6)), # magic
               list(n=as.integer(n)),
               PACKAGE="microsimulation")
  out <- transform(as.data.frame(out),
                   state=enum(state,stateT),
                   event=enum(event,eventT))
  ## tidy up
  state$reset()
  out
}

callSimplePerson <- function(n=10) {
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call("callSimplePerson",
               parms=list(n=as.integer(n)),
               PACKAGE="microsimulation")
  out <- transform(as.data.frame(out),
                   state=enum(state,stateT),
                   event=enum(event,eventT))
  out
}

callSimplePerson2 <- function(n=10) {
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call("callSimplePerson2",
               parms=list(n=as.integer(n)),
               PACKAGE="microsimulation")
  out$events$state <- enum(out$events$state,stateT)
  out$events$event <- enum(out$events$event,eventT)
  out$pt$state <- enum(out$pt$state,stateT)
  out$prev$state <- enum(out$prev$state,stateT)
  out
}


.testPackage <- function() {
  list(callSimplePerson(),
       callPersonSimulation(n=10),
       callSimplePerson2())
}

EventQueue <-
  setRefClass("EventQueue",
              fields = list(times = "numeric", events = "list"),
              methods = list(
                insert = function(time,event) {
                  ord <- order(newtimes <- c(times,time))
                  times <<- newtimes[ord]
                  events <<- c(events,list(event))[ord]
                },
                pop = function() {
                  head <- structure(events[[1]], time=times[1])
                  times <<- times[-1]
                  events <<- events[-1]
                  return(head)
                },
                empty = function() length(times) == 0,
                clear = function() {
                  times <<- numeric()
                  events <<- list()
                },
                remove = function(predicate) {
                  i <- sapply(events, predicate)
                  stopifnot(is.logical(i))
                  i[is.na(i)] <- TRUE
                  times <<- times[!i]
                  events <<- events[!i]
                }))

BaseDiscreteEventSimulation <-
  setRefClass("BaseDiscreteEventSimulation",
              contains = "EventQueue",
              fields = list(currentTime = "numeric",
                previousEventTime = "numeric"),
              methods = list(
                scheduleAt = function(time, event) insert(time, event),
                init = function() stop("VIRTUAL!"), 
                handleMessage = function(event) stop("VIRTUAL!"),
                final = function() {},
                now = function() currentTime,
                run = function() {
                  previousEventTime <<- 0.0
                  init()
                  while(!empty()) {
                    event <- pop()
                    currentTime <<- attr(event, "time")
                    handleMessage(event)
                    previousEventTime <<- currentTime
                  }
                  final()
                }))

RNGStream <- function(nextStream = TRUE, iseed = NULL) {
  stopifnot(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) || !is.null(iseed))
  oldseed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else NULL
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  if (!is.null(iseed)) set.seed(iseed)
  current <- if (nextStream) parallel::nextRNGStream(.Random.seed) else .Random.seed
  .Random.seed <<- startOfStream <- startOfSubStream <- current
  structure(list(resetRNGkind = function() {
    if (!is.null(oldseed)) 
      assign(".Random.seed", oldseed, envir = .GlobalEnv)
    else rm(.Random.seed, envir = .GlobalEnv)
  },
                 open = function() .Random.seed <<- current,
                 close = function() current <<- .Random.seed,
                 resetStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream,
                 resetSubStream = function() .Random.seed <<- current <<- startOfSubStream,
                 nextSubStream = function() .Random.seed <<- current <<- startOfSubStream <<- parallel::nextRNGSubStream(startOfSubStream),
                 nextStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream <<- parallel::nextRNGStream(startOfStream)),
            class="RNGStream")
  }

with.RNGStream <- function(stream,expr,...) {
  stream$open()
  out <- eval(substitute(expr), enclos = parent.frame(), ...)
  stream$close()
  out
}
setOldClass("RNGStream")
