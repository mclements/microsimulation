
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

callPersonSimulation <- function(seed=rep(12345,6),
                                 n=500) {
  stateT =c("Healthy","Localised","DxLocalised","LocallyAdvanced",
    "DxLocallyAdvanced","Metastatic","DxMetastatic","Death")
  eventT =c("toDeath", "toPCDeath", "toLocalised", "toDxLocalised",
	      "toDxLocallyAdvanced",
	      "toLocallyAdvanced", "toMetastatic", "toDxMetastatic")
  out <- .Call("callPersonSimulation",
               as.integer(seed),
               list(n=as.integer(n)),
               PACKAGE="microsimulation")
  out <- transform(as.data.frame(out),
                   state=enum(state,stateT),
                   event=enum(event,eventT))
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

RNGStream <- function(nextStream = TRUE) {
  current <- if (nextStream) nextRNGStream(.Random.seed) else .Random.seed
  .Random.seed <<- startOfStream <- startOfSubStream <- current
  structure(list(open = function() .Random.seed <<- current,
                 close = function() current <<- .Random.seed,
                 resetStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream,
                 resetSubStream = function() .Random.seed <<- current <<- startOfSubStream,
                 nextSubStream = function() .Random.seed <<- current <<- startOfSubStream <<- nextRNGSubStream(startOfSubStream),
                 nextStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream <<- nextRNGStream(startOfStream)),
            class="RNGStream")
  }

with.RNGStream <- function(stream,expr,...) {
  stream$open()
  on.exit(stream$close())
  eval(substitute(expr), enclos = parent.frame(), ...)
}
setOldClass("RNGStream")
