
## initial values for the cervical model
CervicalParameters <- list(
    nLifeHistories = 10L, screen = 0L, ## integers
    discountRate.effectiveness = 0.03,
    discountRate.costs = 0.03,
    full_report = 1.0,
    ## this needs to be changed
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
    cost_parameters = c(Invitation = 50,
        FormalPSA = 130,
        OpportunisticPSA = 1910,
        FormalPSABiomarker = 730,
        OpportunisticPSABiomarker = 2510, #N.B. This one is new and should be used
        Biopsy = 12348,
        Prostatectomy = 117171,
        RadiationTherapy = 117171,
        ActiveSurveillance = 141358,
        CancerDeath = 585054, 
        Death = 0),
    ## IHE doesn't use the postrecovery period (as reported in the Heijnsdijk 2012 reference), should we?
    utility_estimates = 1 - c(Invitation = 1,
        FormalPSA = 0.99,
        FormalPSABiomarker = 0.90,
        Biopsy = 0.90,
        OpportunisticPSA = 0.99,
        ProstatectomyPart1 = 0.67,
        ProstatectomyPart2 = 0.77,
        RadiationTherapyPart1 = 0.73,
        RadiationTherapyPart2 = 0.78,
        ActiveSurveillance = 0.97,
        PalliativeTherapy = 0.60,
        TerminalIllness = 0.40,
        MetastaticCancer = 0.85,
        Death = 0.00),
    ## Utility duration is given in years.
    utility_duration = c(Invitation = 0.0,
        FormalPSA = 1/52,
        FormalPSABiomarker = 3/52,
        Biopsy = 3/52,
        OpportunisticPSA = 1/52,
        ProstatectomyPart1 = 2/12,
        ProstatectomyPart2 = 10/12,
        RadiationTherapyPart1 = 2/12,
        RadiationTherapyPart2 = 10/12,
        ActiveSurveillance = 7,
        PalliativeTherapy = 30/12,
        TerminalIllness = 6/12)
    )
## This needs to be changed
pop1 <- data.frame(cohort=2012:1900,
                   pop=c(rep(17239,9), 16854, 16085, 15504, 15604, 16381, 16705, 
                       16762, 16853, 15487, 14623, 14066, 13568, 13361, 13161, 13234, 
                       13088, 12472, 12142, 12062, 12078, 11426, 12027, 11963, 12435, 
                       12955, 13013, 13125, 13065, 12249, 11103, 9637, 9009, 8828, 
                       8350, 7677, 7444, 7175, 6582, 6573, 6691, 6651, 6641, 6268, 
                       6691, 6511, 6857, 7304, 7308, 7859, 7277, 8323, 8561, 7173, 
                       6942, 7128, 6819, 5037, 6798, rep(6567,46)))
CervicalData <- list()
cervicalEnum <- list(stateT=NULL, eventT=NULL)
## Hcervical <- data.frame(hpv, age, from, to, survival)
Hcervical <- data.frame()
hpvT <- c("LR_HPV","HPV_16","HPV_18","Other_HR_HPV")
cervStateT <- c("Normal", "HPV", "CIN1", "CIN23", "LocalCancer", "RegionalCancer", "DistantCancer", "Death")
cervEventT <- c("toHPV", "toCIN1", "toNormal", "toCIN23", "toNoCIN", "toLocalCancer", "toRegionalCancer,
		toDistantCancer", "toUtility", "toUtilityChange", "toOtherDeath", "toCancerDeath")
cervicalData <- list(H=Hcervical)

callCervical <- function(n=10, nLifeHistories=10,
                      seed=12345,
                      flatPop = FALSE, pop = pop1, tables = list(), debug=FALSE,
                      discountRate = 0.03, parms = NULL, mc.cores=1) {
  ## save the random number state for resetting later
  state <- RNGstate(); on.exit(state$reset())
  ## yes, we use the user-defined RNG
  RNGkind("user")
  set.user.Random.seed(seed)
  ## birth cohorts that should give approximately the number of men alive in Stockholm in 2012
  ## check the input arguments
  stopifnot(is.na(n) || is.integer(as.integer(n)))
  stopifnot(is.integer(as.integer(nLifeHistories)))
  ## NB: sample() calls the random number generator (!)
  if (is.vector(pop)) {
      flatPop <- TRUE
      pop <- data.frame(cohort=pop,pop=1)
  }
  if (is.na(n)) {
    cohort <- pop$cohort[rep.int(1:nrow(pop),times=pop$pop)]
    n <- length(cohort)
  } else {
      if (flatPop) {
          cohort <- rep(pop$cohort,each=ceiling(n/nrow(pop))) #Need ceiling so int n=!0
          n <- ceiling(n/nrow(pop)) * nrow(pop) #To get the chunks right
      } else
          cohort <- sample(pop$cohort,n,prob=pop$pop/sum(pop$pop),replace=TRUE)
  }
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
  ## Minor changes to cervicalData
  if (!is.null(tables)) 
      for (name  in names(tables))
          cervicalData[[name]] <- tables[[name]]
  updateParameters <- c(parms,
                        list(nLifeHistories=as.integer(nLifeHistories),
                             discountRate.costs=discountRate,
                             discountRate.effectiveness=discountRate))
  parameter <- CervicalParameters
  for (name in names(updateParameters))
      parameter[[name]] <- updateParameters[[name]]
  pind <- sapply(parameter,class)=="numeric" & sapply(parameter,length)==1
  bInd <- sapply(parameter,class)=="logical" & sapply(parameter,length)==1
  ## now run the chunks separately
  print(system.time(out <- parallel::mclapply(1:mc.cores,
                function(i) {
                  chunk <- chunks[[i]]
                  set.user.Random.seed(initialSeeds[[i]])
                  .Call("callCervical",
                        parms=list(n=as.integer(length(chunk)),
                            firstId=ns[i],
                            debug=debug, # bool
                            cohort=as.double(chunk),
                            parameter=unlist(parameter[pind]),
                            bparameter=unlist(parameter[bInd]),
                            otherParameters=parameter[!pind & !bInd],
                            tables=cervicalData),
                        PACKAGE="microsimulation")
                }, mc.cores = mc.cores)))
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
  if (length(out[[1]]$summary) == 0) summary <- list()
  else {
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
      enum(summary$events$event) <- eventT
}
  map2df <- function(obj) as.data.frame(do.call("cbind",obj))
  lifeHistories <- do.call("rbind",lapply(out,function(obj) map2df(obj$lifeHistories)))
  parameters <- map2df(out[[1]]$parameters)
  ## Identifying elements without name which also need to be rbind:ed
  ## costs <- do.call("rbind",lapply(out,function(obj) data.frame(obj$costs)))
  ## names(costs) <- c("item","cohort","age","costs")
  names(lifeHistories) <- c("id","state","ext_grade","dx","event","begin","end","year","psa")
  enum(lifeHistories$state) <- stateT
  enum(lifeHistories$event) <- eventT
  enum <- list(stateT = stateT, eventT = eventT)
  out <- list(n=n,enum=enum,lifeHistories=lifeHistories,
              parameters=parameters,
              ## summary=summary, costs=costs,
              cohort=data.frame(table(cohort)),
              discountRate = discountRate)
  class(out) <- "cervical"
  out
}
