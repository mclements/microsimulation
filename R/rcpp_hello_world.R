
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

## initial values for the FHCRC model
FhcrcParameters <- list(
    tau2 = 0.0829,
    g0=0.0005,
    gm=0.0004,
    gc=0.0015, 
    thetac=19.1334,
    mubeta0=-1.609,
    sebeta0=0.2384,
    mubeta1=0.04463,
    sebeta1=0.0430,
    mubeta2=c(0.0397,0.1678),
    sebeta2=c(0.0913,0.3968),
    c_txlt_interaction = 1.0,
    c_baseline_specific = 1.0,
    screeningCompliance = 0.75,
    biopsyCompliance = 0.858,
    biopsySensitivity = 0.8,
    studyParticipation = 35.0/260.0,
    nLifeHistories = 10L, screen = 0L, ## integers
    psaThreshold = 3.0,
    psaThresholdBiopsyFollowUp = 4.0,
    c_low_grade_slope=-0.006,
    discountRate = 0.035,
    mu0=c(0.00219, 0.000304, 5.2e-05, 0.000139, 0.000141, 3.6e-05, 7.3e-05, 
        0.000129, 3.8e-05, 0.000137, 6e-05, 8.1e-05, 6.1e-05, 0.00012, 
        0.000117, 0.000183, 0.000185, 0.000397, 0.000394, 0.000585, 0.000448, 
        0.000696, 0.000611, 0.000708, 0.000659, 0.000643, 0.000654, 0.000651, 
        0.000687, 0.000637, 0.00063, 0.000892, 0.000543, 0.00058, 0.00077, 
        0.000702, 0.000768, 0.000664, 0.000787, 0.00081, 0.000991, 9e-04, 
        0.000933, 0.001229, 0.001633, 0.001396, 0.001673, 0.001926, 0.002217, 
        0.002562, 0.002648, 0.002949, 0.002729, 0.003415, 0.003694, 0.004491, 
        0.00506, 0.004568, 0.006163, 0.006988, 0.006744, 0.00765, 0.007914, 
        0.009153, 0.010231, 0.011971, 0.013092, 0.013839, 0.015995, 0.017693, 
        0.018548, 0.020708, 0.022404, 0.02572, 0.028039, 0.031564, 0.038182, 
        0.042057, 0.047361, 0.05315, 0.058238, 0.062619, 0.074934, 0.089776, 
        0.099887, 0.112347, 0.125351, 0.143077, 0.153189, 0.179702, 0.198436, 
        0.240339, 0.256215, 0.275103, 0.314157, 0.345252, 0.359275, 0.41768, 
        0.430279, 0.463636, 0.491275, 0.549738, 0.354545, 0.553846, 0.461538, 
        0.782609),
    cost_parameters = c(InvitationCost = 15,
        FormalPSACost = 41,
        FormalPSABiomarkerCost = 641,
        BiopsyCost = 8082,
        OpportunisticPSACost = 1774,
        ProstatectomyCost = 95000,
        RadiationTherapyCost = 135000,
        ActiveSurveillanceCost = 140000,
        MetastaticCancerCost = 769574,
        DeathCost = 0),
    ## IHE doesn't use the postrecovery period (as reported in the Heijnsdijk 2012 reference), should we?
    utility_estimates = 1 - c(InvitationUtility = 1,
        FormalPSAUtility = 0.99,
        FormalPSABiomarkerUtility = 0.90,
        BiopsyUtility = 0.90,
        OpportunisticPSAUtility = 0.99,
        ProstatectomyUtilityPart1 = 0.67,
        ProstatectomyUtilityPart2 = 0.77,
        RadiationTherapyUtilityPart1 = 0.73,
        RadiationTherapyUtilityPart2 = 0.78,
        ActiveSurveillanceUtility = 0.97,
        MetastaticCancerUtilityPart1 = 0.60,
        MetastaticCancerUtilityPart2 = 0.40,
        DeathUtility = 0.00),
    ## Utility duration is given in years.
    utility_duration = c(InvitationUtilityDuration = 0.0,
        FormalPSAUtilityDuration = 1/52,
        FormalPSABiomarkerUtilityDuration = 3/52,
        BiopsyUtilityDuration = 3/52,
        OpportunisticPSAUtilityDuration = 1/52,
        ProstatectomyUtilityDurationPart1 = 2/12,
        ProstatectomyUtilityDurationPart2 = 10/12,
        RadiationTherapyUtilityDurationPart1 = 2/12,
        RadiationTherapyUtilityDurationPart2 = 10/12,
        ActiveSurveillanceUtilityDuration = 7,
        MetastaticCancerUtilityDurationPart1 = 30/12,
        MetastaticCancerUtilityDurationPart2 = 6/12)
    )
ParameterNV <- FhcrcParameters[sapply(FhcrcParameters,class)=="numeric" & sapply(FhcrcParameters,length)==1]
## ParameterIV <- FhcrcParameters[sapply(FhcrcParameters,class)=="integer" & sapply(FhcrcParameters,length)==1]


callFhcrc <- function(n=10,screen="noScreening",nLifeHistories=10,screeningCompliance=0.75,
                      seed=12345, studyParticipation=50/260, psaThreshold=3.0, includePSArecords=FALSE, mc.cores=1) {
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
              "toCancerDeath","toOtherDeath","toScreen",
              "toScreenInitiatedBiopsy","toClinicalDiagnosticBiopsy","toScreenDiagnosis",
              "toOrganised","toTreatment","toCM","toRP","toRT","toADT","toChangeUtility","toAgeUtility")
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
  updateParameters <- list(nLifeHistories=as.integer(nLifeHistories),
                           screeningCompliance=screeningCompliance,
                           studyParticipation=studyParticipation,
                           psaThreshold=psaThreshold,
                           screen=as.integer(screenIndex))
  parameter <- FhcrcParameters
  for (name in names(updateParameters))
      parameter[[name]] <- updateParameters[[name]]
  pind <- sapply(parameter,class)=="numeric" & sapply(parameter,length)==1
  ## now run the chunks separately
  print(system.time(out <- parallel::mclapply(1:mc.cores,
                function(i) {
                  chunk <- chunks[[i]]
                  set.user.Random.seed(initialSeeds[[i]])
                  .Call("callFhcrc",
                        parms=list(n=as.integer(length(chunk)),
                            firstId=ns[i],
                            cohort=as.double(chunk),
                            parameter=unlist(parameter[pind]),
                            otherParameters=parameter[!pind],
                            tables=fhcrcData,
                            includePSArecords=includePSArecords),
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
  rbindList <- function(obj) # recursive  
      if (is.list(obj)) do.call("rbind",lapply(obj,rbindList)) else data.frame(obj)
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
  ## grab all of the pt, prev, ut, events from summary
  ## pt <- lapply(out, function(obj) obj$summary$pt)
  summary <- lapply(seq_along(out[[1]]$summary),
                    function(i) do.call("rbind",
                                        lapply(out, function(obj) reader(obj$summary[[i]]))))
  names(summary) <- names(out[[1]]$summary)
  states <- c("state","grade","dx","psa","cohort")
  names(summary$prev) <- c(states,"age","count")
  names(summary$pt) <- c(states,"age","pt")
  names(summary$ut) <- c(states,"age","ut")
  names(summary$events) <- c(states,"event","age","n")
  summary <- lapply(summary,function(obj) within(obj,year <- cohort+age))
  ## map2df <- function(obj) "names<-"(data.frame(obj[-1]),obj[[1]]) 
  map2df <- function(obj) as.data.frame(do.call("cbind",obj))
  lifeHistories <- do.call("rbind",lapply(out,function(obj) map2df(obj$lifeHistories)))
  psarecord <- do.call("rbind",lapply(out,function(obj) data.frame(obj$psarecord)))
  parameters <- map2df(out[[1]]$parameters)
  ## Identifying elements without name which also need to be rbind:ed
  costsNameless_idx <- names(out[[1]]$costs)==""
  costs <- cbind(rbindList(out[[1]]$costs[costsNameless_idx]),cbindList(out[[1]]$costs[!costsNameless_idx]))
  names(costs) <- c("item","age","costs")
  names(lifeHistories) <- c("id","state","ext_grade","dx","event","begin","end","psa")
  enum(summary$events$event) <- eventT
  enum(lifeHistories$state) <- stateT
  enum(lifeHistories$dx) <- diagnosisT
  enum(lifeHistories$event) <- eventT
  enum <- list(stateT = stateT, eventT = eventT, screenT = screenT, diagnosisT = diagnosisT,
               psaT = psaT)
  out <- list(n=n,screen=screen,enum=enum,lifeHistories=lifeHistories,parameters=parameters,
              ## prev=summary$prev, pt=summary$pt, events=summary$events)
              summary=summary,costs=costs, psarecord=psarecord)
  class(out) <- "fhcrc"
  out
}

## R --slave -e "options(width=200); require(microsimulation); callFhcrc(100,nLifeHistories=1e5,screen=\"screen50\")[[\"parameters\"]]"


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
