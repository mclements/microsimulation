
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

set.user.Random.seed <- function (seed) {
  seed <- as.integer(seed)
  stopifnot(is.integer(seed))
  if (length(seed) == 1) seed <- rep(seed,6)
  .C("r_set_user_random_seed",seed = seed,PACKAGE="microsimulation")
  return(invisible(TRUE))
}

user.Random.seed <- function() {
  .C("r_get_user_random_seed", seed=rep(1L,6), PACKAGE="microsimulation")
}

enum <- function(obj, labels)
  factor(obj, levels=0:(length(labels)-1), labels=labels)

"enum<-" <- function(obj, value) {
  enum(if(is.factor(obj)) unclass(obj)-1 else obj, value)
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
  }
  list(oldseed = oldseed, reset = reset)
}
  
callPersonSimulation <- function(n=20,seed=rep(12345,6)) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("user")
  set.user.Random.seed(seed)
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
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath", "toCheck")
  out <- .Call("callSimplePerson",
               parms=list(n=as.integer(n)),
               PACKAGE="microsimulation")
  out <- transform(as.data.frame(out),
                   state=enum(state,stateT),
                   event=enum(event,eventT))
  out
}

callSimplePerson2 <- function(n=10) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call("callSimplePerson2",
               parms=list(n=as.integer(n)),
               PACKAGE="microsimulation")
  reader <- function(obj) cbind(data.frame(state=enum(obj$state[[1]],stateT)),data.frame(obj[-1]))
  out <- lapply(out,reader)
  enum(out$events$event) <- eventT
  out
}

callFhcrc <- function(n=10,screen="noScreening",nLifeHistories=10,screeningCompliance=0.75,seed=12345) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("user")
  set.user.Random.seed(seed)
  ## birth cohorts that should give approximately the number of men alive in Stockholm in 2012
  pop1 <- data.frame(cohort=1972:1900, pop=c(17239, 16854, 16085, 15504, 15604, 16381, 16705, 
    16762, 16853, 15487, 14623, 14066, 13568, 13361, 13161, 13234, 
    13088, 12472, 12142, 12062, 12078, 11426, 12027, 11963, 12435, 
    12955, 13013, 13125, 13065, 12249, 11103, 9637, 9009, 8828, 
    8350, 7677, 7444, 7175, 6582, 6573, 6691, 6651, 6641, 6268, 
    6691, 6511, 6857, 7304, 7308, 7859, 7277, 8323, 8561, 7173, 
    6942, 7128, 6819, 5037, 6798, rep(6567,14)))
  cohort <- pop1$cohort[rep.int(1:nrow(pop1),times=pop1$pop)]
  screenT <- c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70", "fourYearlyScreen50to70", "screen50",
               "screen60", "screen70", "screenUptake", "stockholm3_goteborg",
               "stockholm3_risk_stratified")
  stateT <- c("Healthy","Localised","Metastatic")
  eventT <- c("toLocalised","toMetastatic","toClinicalDiagnosis",
              "toCancerDeath","toOtherDeath","toScreen","toBiopsy","toScreenDiagnosis",
              "toOrganised")
  diagnosisT <- c("NotDiagnosed","ClinicalDiagnosis","ScreenDiagnosis")
  psaT <- c("PSA<3","PSA>=3")
  stopifnot(screen %in% screenT)
  stopifnot(is.na(n) || is.integer(as.integer(n)))
  stopifnot(is.integer(as.integer(nLifeHistories)))
  stopifnot(is.double(as.double(screeningCompliance)))
  if (is.na(n)) n <- length(cohort) else cohort <- sample(pop1$cohort,n,prob=pop1$pop/sum(pop1$pop),replace=TRUE)
  screenIndex <- which(screen == screenT) - 1
  out <- .Call("callFhcrc",
               parms=list(n=as.integer(n),screen=as.integer(screenIndex),nLifeHistories=as.integer(nLifeHistories),
                 screeningCompliance=as.double(screeningCompliance),
                 cohort=as.double(cohort)),
               PACKAGE="microsimulation")
  reader <- function(obj) {
    out <- cbind(data.frame(state=enum(obj$state[[1]],stateT),
                            dx=enum(obj$state[[2]],diagnosisT),
                            psa=enum(obj$state[[3]],psaT),
                            cohort=obj$state[[4]]),
                 data.frame(obj[-1]))
    out$year <- out$cohort + out$age
    out
  }
  out$summary <- lapply(out$summary,reader)
  enum(out$summary$events$event) <- eventT
  enum(out$lifeHistories$state) <- stateT
  enum(out$lifeHistories$dx) <- diagnosisT
  enum(out$lifeHistories$event) <- eventT
  out$lifeHistories <- data.frame(out$lifeHistories)
  out$parameters <- data.frame(out$parameters)
  out$enum <- list(stateT = stateT, eventT = eventT, screenT = screenT, diagnosisT = diagnosisT,
                   psaT = psaT)
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
