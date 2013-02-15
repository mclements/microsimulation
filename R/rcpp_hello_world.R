
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

"enum<-" <- function(obj, value) {
  enum(obj,value)
}

RNGstate <- function() {
  ## house-keeping for random streams (see parallel::clusterSetRNGStream)
  oldseed <- if (exists(".Random.seed", envir = .GlobalEnv, 
                        inherits = FALSE)) 
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else NULL
  reset <- function() {
    if (!is.null(oldseed)) 
      assign(".Random.seed", oldseed, envir = .GlobalEnv)
    else rm(.Random.seed, envir = .GlobalEnv)
  }
  list(oldseed = oldseed, reset = reset)
}
  
callPersonSimulation <- function(n=20,seed=rep(12345,6)) {
  state <- RNGstate(); on.exit(state$reset())
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
  enum(out$events$state) <- stateT
  enum(out$events$event) <- eventT
  enum(out$pt$state) <- stateT
  enum(out$prev$state) <- stateT
  out
}

callFhcrc <- function(n=10,screen="noScreening",nLifeHistories=10,screeningCompliance=0.5) {
  screenT <- c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70", "fourYearlyScreen50to70", "screen50",
               "screen60", "screen70")
  stateT <- c("Healthy","Localised","Metastatic","ClinicalDiagnosis","ClinicalMetastaticDiagnosis","ScreenDiagnosis","ScreenMetastaticDiagnosis")
  eventT <- c("toLocalised","toMetastatic","toClinicalDiagnosis","toClinicalMetastaticDiagnosis",
              "toCancerDeath","toOtherDeath","toScreen","toBiopsy","toScreenDiagnosis","toScreenMetastaticDiagnosis")
  stopifnot(screen %in% screenT)
  screenIndex <- which(screen == screenT) - 1
  out <- .Call("callFhcrc",
               parms=list(n=as.integer(n),screen=as.integer(screenIndex),nLifeHistories=as.integer(nLifeHistories),
                 screeningCompliance=as.double(screeningCompliance)),
               PACKAGE="microsimulation")
  enum(out$summary$events$state) <- stateT
  enum(out$summary$events$event) <- eventT
  enum(out$summary$pt$state) <- stateT
  enum(out$summary$prev$state) <- stateT
  enum(out$lifeHistories$state) <- stateT
  enum(out$lifeHistories$event) <- eventT
  out$lifeHistories <- data.frame(out$lifeHistories)
  out$parameters <- data.frame(out$parameters)
  out$enum <- list(stateT = stateT, eventT = eventT, screenT = screenT)
  out$n <- n
  out
}

callFhcrcTest <- function(n=10,screen="noScreening",nLifeHistories=10,screeningCompliance=0.5) {
  screenT <- c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70", "fourYearlyScreen50to70", "screen50",
               "screen60", "screen70")
  stateT <- c("Healthy","Localised","Metastatic")
  eventT <- c("toLocalised","toMetastatic","toClinicalDiagnosis",
              "toCancerDeath","toOtherDeath","toScreen","toBiopsy","toScreenDiagnosis")
  diagnosisT <- c("NotDiagnosed","ClinicalDiagnosis","ScreenDiagnosis")
  stopifnot(screen %in% screenT)
  screenIndex <- which(screen == screenT) - 1
  out <- .Call("callFhcrcTest",
               parms=list(n=as.integer(n),screen=as.integer(screenIndex),nLifeHistories=as.integer(nLifeHistories),
                 screeningCompliance=as.double(screeningCompliance)),
               PACKAGE="microsimulation")
  enum(out$summary$events$state1) <- stateT
  enum(out$summary$events$state2) <- diagnosisT
  enum(out$summary$events$event) <- eventT
  enum(out$summary$pt$state1) <- stateT
  enum(out$summary$pt$state2) <- diagnosisT
  enum(out$summary$prev$state1) <- stateT
  enum(out$summary$prev$state2) <- diagnosisT
  enum(out$lifeHistories$state) <- stateT
  enum(out$lifeHistories$dx) <- diagnosisT
  enum(out$lifeHistories$event) <- eventT
  out$lifeHistories <- data.frame(out$lifeHistories)
  out$parameters <- data.frame(out$parameters)
  out$enum <- list(stateT = stateT, eventT = eventT, screenT = screenT, diagnosisT = diagnosisT)
  out$n <- n
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
