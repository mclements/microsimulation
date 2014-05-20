
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
## http://r.789695.n4.nabble.com/How-to-construct-a-valid-seed-for-l-Ecuyer-s-method-with-given-Random-seed-td4656340.html
unsigned <- function(seed) ifelse(seed < 0, seed + 2^32, seed)
signed <- function(seed) ifelse(seed>2^31, seed-2^32, seed)

set.user.Random.seed <- function (seed) {
  seed <- as.double(unsigned(seed))
  if (length(seed) == 1) seed <- rep(seed,6)
  if (length(seed) == 7) seed <- seed[-1]
  .C("r_set_user_random_seed",seed = seed,PACKAGE="microsimulation")
  return(invisible(seed))
}

next.user.Random.substream <- function () {
  .C("r_next_rng_substream",PACKAGE="microsimulation")
  return(invisible(TRUE))
}

user.Random.seed <- function() {
  c(407L,
    as.integer(signed(.C("r_get_user_random_seed", seed=as.double(rep(1,6)),
                         PACKAGE="microsimulation")$seed)))
}

enum <- function(obj, labels, start=0) {
  if (is.logical(obj)) obj <- obj+0
  structure(factor(obj, levels=start + (0:(length(labels)-1)), labels=labels),start=start)
}

"enum<-" <- function(obj, value) {
  enum(if(is.factor(obj)) unclass(obj)-1+attr(obj,"start") else obj, value)
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
    else {
        ## clean up if we created a .Random.seed
        if (exists(".Random.seed" ,envir = .GlobalEnv, inherits = FALSE))
            rm(.Random.seed, envir = .GlobalEnv, inherits = FALSE)
    }
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
  out <- `names<-`(as.data.frame(out$Value),out$Key)
  out <- transform(out,
                   state=enum(state,stateT),
                   event=enum(event,eventT))
  ## tidy up
  out
}

callSimplePerson <- function(n=10) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  set.seed(12345)
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
  set.seed(12345)
  stateT <- c("Healthy","Cancer","Death")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call("callSimplePerson2",
               parms=list(n=as.integer(n)),
               PACKAGE="microsimulation")
  reader <- function(obj) 
      cbind(data.frame(state=enum(obj[[1]],stateT)),
          data.frame(obj[-1]))
  out <- lapply(out,reader)
  out$events <- with(out$events, data.frame(state=state,event=enum(event,eventT),age=age,number=number))
  out$pt <- with(out$pt, data.frame(state=state,age=age,pt=pt))
  out$prev <- with(out$prev, data.frame(state=state,age=age,number=number))
  out
}


## readEventReport <- function(obj) {
##     pt <- obj$pt
##     n <- ncol(pt)
##     names(pt)[(n-1):n] <- c("age","pt")
## }

callIllnessDeath <- function(n=10L,cure=0.1,zsd=0) {
  state <- RNGstate(); on.exit(state$reset())
  RNGkind("Mersenne-Twister")
  set.seed(12345)
  stateT <- c("Healthy","Cancer")
  eventT <- c("toOtherDeath", "toCancer", "toCancerDeath")
  out <- .Call("callIllnessDeath",
               parms=list(n=as.integer(n),cure=as.double(cure),zsd=as.double(zsd)),
               PACKAGE="microsimulation")
  reader <- function(obj) 
      cbind(data.frame(state=enum(obj[[1]],stateT)),
          data.frame(obj[-1]))
  out <- lapply(out,reader)
  out$events <- with(out$events, data.frame(state=state,event=enum(event,eventT),age=age,number=number))
  out$pt <- with(out$pt, data.frame(state=state,age=age,pt=pt))
  out$prev <- with(out$prev, data.frame(state=state,age=age,number=number))
  out
}

callFhcrc <- function(n=10,screen="noScreening",nLifeHistories=10,screeningCompliance=0.75,
                      seed=12345, studyParticipation=50/260, psaThreshold=3.0, mc.cores=1) {
  ## save the random number state for resetting later
  state <- RNGstate(); on.exit(state$reset())
  ## yes, we use the user-defined RNG
  RNGkind("user")
  set.user.Random.seed(seed)
  ## birth cohorts that should give approximately the number of men alive in Stockholm in 2012
  pop1 <- data.frame(cohort=1980:1900, pop=c(rep(17239,9), 16854, 16085, 15504, 15604, 16381, 16705, 
    16762, 16853, 15487, 14623, 14066, 13568, 13361, 13161, 13234, 
    13088, 12472, 12142, 12062, 12078, 11426, 12027, 11963, 12435, 
    12955, 13013, 13125, 13065, 12249, 11103, 9637, 9009, 8828, 
    8350, 7677, 7444, 7175, 6582, 6573, 6691, 6651, 6641, 6268, 
    6691, 6511, 6857, 7304, 7308, 7859, 7277, 8323, 8561, 7173, 
    6942, 7128, 6819, 5037, 6798, rep(6567,14)))
  ## these enum strings should be moved to C++
  screenT <- c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70", "fourYearlyScreen50to70", "screen50",
               "screen60", "screen70", "screenUptake", "stockholm3_goteborg",
               "stockholm3_risk_stratified")
  stateT <- c("Healthy","Localised","Metastatic")
  gradeT <- c("Gleason_le_6","Gleason_7","Gleason_ge_8")
  eventT <- c("toLocalised","toMetastatic","toClinicalDiagnosis",
              "toCancerDeath","toOtherDeath","toScreen","toBiopsy","toScreenDiagnosis",
              "toOrganised","toTreatment","toCM","toRP","toRT","toADT")
  diagnosisT <- c("NotDiagnosed","ClinicalDiagnosis","ScreenDiagnosis")
  treatmentT <- c("CM","RP","RT")
  psaT <- c("PSA<3","PSA>=3") # not sure where to put this...
  ## check the input arguments
  stopifnot(screen %in% screenT)
  stopifnot(is.na(n) || is.integer(as.integer(n)))
  stopifnot(is.integer(as.integer(nLifeHistories)))
  stopifnot(is.double(as.double(screeningCompliance)))
  screenIndex <- which(screen == screenT) - 1
  ## NB: sample() calls the random number generator (!)
  if (is.na(n)) {
    cohort <- pop1$cohort[rep.int(1:nrow(pop1),times=pop1$pop)]
    n <- length(cohort)
  } else
    cohort <- sample(pop1$cohort,n,prob=pop1$pop/sum(pop1$pop),replace=TRUE)
  cohort <- sort(cohort)
  ## now separate the data into chunks
  chunks <- tapply(cohort, sort((0:(n-1)) %% mc.cores), I)
  ## set the initial random numbers
  currentSeed <- user.Random.seed()
  powerFun <- function(obj,FUN,n,...) {
    for(i in 1:n)
      obj <- FUN(obj,...)
    obj
  }
  initialSeeds <- Reduce(function(seed,i) powerFun(seed,parallel::nextRNGStream,10),
                         1:mc.cores, currentSeed, accumulate=TRUE)[-1]
  ns <- cumsum(sapply(chunks,length))
  ns <- c(0,ns[-length(ns)])
  ## Minor changes to fhcrcData
  fhcrcData$prtx$Age <- as.double(fhcrcData$prtx$Age)
  fhcrcData$prtx$DxY <- as.double(fhcrcData$prtx$DxY)
  fhcrcData$prtx$G <- fhcrcData$prtx$G - 1L
  fhcrcData$pradt$Grade <- fhcrcData$pradt$Grade - 1L
  fhcrcData$pradt$Age <- as.double(fhcrcData$pradt$Age)
  fhcrcData$pradt$DxY <- as.double(fhcrcData$pradt$DxY)
  fhcrcData$survival_local <-
      with(fhcrcData$survival_local,
           data.frame(Age=as.double(AgeLow),Grade=Grade,Time=as.double(Time),
                      Survival=Survival))
  fhcrcData$survival_dist <-
      with(fhcrcData$survival_dist,
           data.frame(Grade=Grade,Time=as.double(Time),
                      Survival=Survival))
  ## now run the chunks separately
  print(system.time(out <- parallel::mclapply(1:mc.cores,
                function(i) {
                  chunk <- chunks[[i]]
                  set.user.Random.seed(initialSeeds[[i]])
                  .Call("callFhcrc",
                        parms=list(n=as.integer(length(chunk)),
                            firstId=ns[i],
                          screen=as.integer(screenIndex),
                          nLifeHistories=as.integer(nLifeHistories),
                          screeningCompliance=as.double(screeningCompliance),
                          studyParticipation=as.double(studyParticipation),
                          psaThreshold=as.double(psaThreshold),
                          cohort=as.double(chunk),
                          tables=fhcrcData),
                        PACKAGE="microsimulation")
                })))
  ## Apologies: we now need to massage the chunks from C++
  ## reader <- function(obj) {
  ##   out <- cbind(data.frame(state=enum(obj$state[[1]],stateT),
  ##                           dx=enum(obj$state[[2]],diagnosisT),
  ##                           psa=enum(obj$state[[3]],psaT),
  ##                           cohort=obj$state[[4]]),
  ##                data.frame(obj[-1]))
  ##   out$year <- out$cohort + out$age
  ##   out
  ## }
  cbindList <- function(obj) # recursive
    if (is.list(obj)) do.call("cbind",lapply(obj,cbindList)) else data.frame(obj)
  reader <- function(obj) {
    obj <- cbindList(obj)
    out <- cbind(data.frame(state=enum(obj[[1]],stateT),
                            grade=enum(obj[[2]],gradeT),
                            dx=enum(obj[[3]],diagnosisT),
                            psa=enum(obj[[4]],psaT),
                            cohort=obj[[5]]),
                 data.frame(obj[,-(1:5)]))
    out
  }
  summary <- lapply(seq_along(out[[1]]$summary),
                    function(i) do.call("rbind",
                                        lapply(out, function(obj) reader(obj$summary[[i]]))))
  names(summary) <- names(out[[1]]$summary)
  states <- c("state","grade","dx","psa","cohort")
  names(summary$prev) <- c(states,"age","count")
  names(summary$pt) <- c(states,"age","pt")
  names(summary$events) <- c(states,"event","age","n")
  summary <- lapply(summary,function(obj) within(obj,year <- cohort+age))
  ## map2df <- function(obj) "names<-"(data.frame(obj[-1]),obj[[1]]) 
  map2df <- function(obj) as.data.frame(do.call("cbind",obj))
  lifeHistories <- do.call("rbind",lapply(out,function(obj) map2df(obj$lifeHistories)))
  parameters <- map2df(out[[1]]$parameters)
  enum(summary$events$event) <- eventT
  enum(lifeHistories$state) <- stateT
  enum(lifeHistories$dx) <- diagnosisT
  enum(lifeHistories$event) <- eventT
  enum <- list(stateT = stateT, eventT = eventT, screenT = screenT, diagnosisT = diagnosisT,
               psaT = psaT)
  out <- list(n=n,screen=screen,enum=enum,lifeHistories=lifeHistories,parameters=parameters,
              ## prev=summary$prev, pt=summary$pt, events=summary$events)
              summary=summary)
  class(out) <- "fhcrc"
  out
}

print.fhcrc <- function(obj,...)
    cat(sprintf("FHCRC prostate cancer model with %i individual(s) under scenario '%s'.\n",
                obj$n, obj$screen),
        ...)

## utility - not exported
assignList <- function(lst,...)
  for(i in 1:length(lst))
    assign(names(lst)[i], lst[[i]], ...)
## assignList(formals(callFhcrc),pos=1)



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
                remove = function(predicate, ...) {
                  i <- sapply(events, predicate, ...)
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
                  scheduleAt = function(time, event) {
                      attr(event,"time") <- time
                      attr(event,"sendingTime") <- currentTime
                      insert(time, event)
                  },
                init = function() stop("VIRTUAL!"), 
                handleMessage = function(event) stop("VIRTUAL!"),
                final = function() {},
                now = function() currentTime,
                reset = function(startTime = 0.0) {
                    clear()
                    previousEventTime <<- currentTime <<- startTime
                },
                run = function(startTime = 0.0) {
                  reset(startTime)
                  init()
                  while(!empty()) {
                    event <- pop()
                    currentTime <<- attr(event, "time")
                    handleMessage(event)
                    previousEventTime <<- currentTime
                  }
                  final()
                }))

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
            if (event %in% c("Death due to other causes", "Cancer death")) {
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
            if (event %in% c("Death due to other causes", "Cancer death")) {
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
            if (event %in% c("Death due to other causes", "Cancer death")) {
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
            if (event %in% c("Death due to other causes", "Cancer death")) {
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
                 seed = function() current,
                 open = function() .Random.seed <<- current,
                 close = function() current <<- .Random.seed,
                 resetStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream,
                 resetSubStream = function() .Random.seed <<- current <<- startOfSubStream,
                 nextSubStream = function() .Random.seed <<- current <<- startOfSubStream <<- parallel::nextRNGSubStream(startOfSubStream),
                 nextStream = function() .Random.seed <<- current <<- startOfSubStream <<- startOfStream <<- parallel::nextRNGStream(startOfStream)),
            class="RNGStream")
  }

with.RNGStream <- function(data,expr,...) {
  data$open()
  out <- eval(substitute(expr), enclos = parent.frame(), ...)
  data$close()
  out
}
setOldClass("RNGStream")
