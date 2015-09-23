
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

rnormPos <- function(n,mean=0,sd=1,lbound=0) {
    if (length(mean)<n) mean <- rep(mean,length=n)
    if (length(sd)<n) sd <- rep(sd,length=n)
    x <- rnorm(n,mean,sd)
    while(any(i <- which(x<lbound)))
        x[i] <- rnorm(length(i),mean[i],sd[i])
    x
}

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

## find the concave frontier
frontier<-function(x,y)
  {
    ichull <- grDevices::chull(cbind(x,y)) # convex hull
    ichull <- ichull[order(x[ichull])]     # order by x
    xi <- x[ichull]
    yi <- y[ichull]          # subset to convex hull
    include <- sapply(1:length(ichull),
                      function(i)       # establish the frontier
                      all(yi[i] >= yi[ xi<xi[i] ]))
    ichull[include]
  }
lines.frontier <- function(x,y,pch=19,type="b",...) {
    index <- frontier(x,y)
    lines(x[index],y[index],pch=pch,type=type,...)
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
    ## panel=FALSE,
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
    ## mubeta2.scale=1.0, # cf. 2.1
    ## beta.rho=0.62,
    c_txlt_interaction = 1.0,
    c_baseline_specific = 1.0,
    c_benefit_value0 = 10, # (value -> reduction): 0.04 -> 10%; 0.18 -> 20%; 10 -> 28%
    sxbenefit = 1.0,
    c_benefit_type = 1, # 0=stage-shift (=> c_benefit_value0=10), 1=lead-time based (=> c_benefit_value1=0.1)
    c_benefit_value1 = 0.1,
    screeningCompliance = 0.75,
    rescreeningCompliance = 0.95,
    biopsyCompliance = 0.858,
    biopsySensitivity = 0.8,
    studyParticipation = 35.0/260.0,
    nLifeHistories = 10L, screen = 0L, ## integers
    psaThreshold = 3.0,
    psaThresholdBiopsyFollowUp = 4.0,
    ## BPThreshold=3.47,
    ## BPThresholdBiopsyFollowUp=3.47,
    ## BPThreshold=4.69,
    ## BPThresholdBiopsyFollowUp=4.69,
    PSA_FP_threshold_nCa=4.4, # reduce FP in no cancers with PSA threshold
    PSA_FP_threshold_GG6=3.6, # reduce FP in GG 6 with PSA threshold
    BPThreshold=4.2, 
    BPThresholdBiopsyFollowUp=4.2, 
    rTPF=1.0,
    rFPF=0.6,
    c_low_grade_slope=-0.006,
    stockholmTreatment = TRUE,
    discountRate.effectiveness = 0.03,
    discountRate.costs = 0.03,
    full_report = 1.0,
    formal_costs = 1.0,
    formal_compliance = 1.0,
    start_screening = 50.0, # start of organised screening
    stop_screening = 70.0,  # end of organised screening
    screening_interval = 2.0, # screening interval for regular_screening
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
IHE <- list(prtx=data.frame(Age=50.0,DxY=1973.0,G=1:2,CM=0.6,RP=0.26,RT=0.14)) ## assumed constant across ages and periods
ParameterNV <- FhcrcParameters[sapply(FhcrcParameters,class)=="numeric" & sapply(FhcrcParameters,length)==1]
## ParameterIV <- FhcrcParameters[sapply(FhcrcParameters,class)=="integer" & sapply(FhcrcParameters,length)==1]
swedenOpportunisticBiopsyCompliance <- cbind(expand.grid(psa=c(4,10),age=c(50,60,70)),
                                compliance=c(0.7, 0.75, 0.6, 0.7, 0.4, 0.5))
swedenFormalBiopsyCompliance <- cbind(expand.grid(psa=c(4,10),age=c(50,60,70)),
                                compliance=0.9)
stockholmTreatment <-
    data.frame(DxY=2008,
               Age=c(50,50,50,55,55,55,60,60,60,65,65,65,70,70,70,75,75,75,80,80,80,85,85,85),
               G=as.integer(c(6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8,6,7,8)-6),
               CM=c(0.231023,0.044872,0.25,0.333333,0.09542,0.328358,0.409439,0.12825,0.348101,0.479167,0.178182,0.401639,0.689013,0.359143,0.56087,0.876543,0.744444,0.809524,1,0.970711,0.952096,0.9375,1,1),
               RP=c(0.700623,0.815839,0.5,0.553592,0.740111,0.326226,0.49981,0.6866,0.291437,0.409879,0.552483,0.279235,0.210949,0.318374,0.179363,0.041152,0.058652,0.049689,0,0.009763,0.023952,0.0625,0,0),
               RT=c(0.068354,0.13929,0.25,0.113074,0.164469,0.345416,0.090751,0.185151,0.360462,0.110954,0.269335,0.319126,0.100038,0.322482,0.259767,0.082305,0.196903,0.140787,0,0.019526,0.023952,0,0,0))
rescreening <- data.frame(age5 = c(30, 30, 30, 30, 35, 35, 35, 35, 40, 40, 
40, 40, 45, 45, 45, 45, 50, 50, 50, 50, 55, 55, 55, 55, 60, 60, 
60, 60, 65, 65, 65, 65, 70, 70, 70, 70, 75, 75, 75, 75, 80, 80, 
80, 80, 85, 85, 85, 85, 90, 90, 90, 90), total_cat = c(0, 1, 
3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 
1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10, 
0, 1, 3, 10, 0, 1, 3, 10, 0, 1, 3, 10), shape = c(1.04293463521838, 
0.825156738737092, 0.899003822110993, 0.589358680965725, 1.05948319345651, 
0.895913335495644, 0.84619233996788, 0.636706369588122, 1.16028493785731, 
0.97919993225436, 0.68848363112076, 0.667566641635962, 1.25727757690268, 
1.11021543939268, 0.717059585208267, 0.724574127291869, 1.3247305288869, 
1.18510576695807, 0.858171614262325, 0.749998577524088, 1.3340371479645, 
1.20792880820659, 0.920717960779785, 0.852052206426261, 1.32268361303857, 
1.24450960892428, 0.988946105418831, 0.919376727690091, 1.31156083562523, 
1.26972967056357, 1.03454490825271, 0.954911145157684, 1.27126976718353, 
1.26045834296745, 1.06684366805728, 0.968088650561122, 1.18494141305856, 
1.20593215449012, 1.05998296831204, 0.989315997334925, 1.12395322065035, 
1.13733607943079, 1.04483735235382, 0.988385269358649, 1.1220716456468, 
1.09128664878198, 1.01326552242001, 0.948730052972523, 1.08315738475919, 
1.09407460180658, 1.05197207574662, 0.978397483932468), cure = c(0.585284664126963, 
0.545117838495959, 0.402804166722055, 0.14283955595853, 0.398430641565923, 
0.327274451509188, 0.148468826128857, 0.00333583482942569, 0.257384245965907, 
0.208165119313289, 0.157385208407254, 0.0696401481749125, 0.166635358789789, 
0.133819199845029, 0.0610524867384188, 0.0535645284203467, 0.134003946557746, 
0.114600739520109, 0.0537234127141967, 0.0308512461181499, 0.104553332825899, 
0.0862550079726627, 0.0399447271614142, 0.0259518317972345, 0.10432761596883, 
0.0823122881511086, 0.0355479853523502, 0.0190476453614591, 0.103738700348739, 
0.0816458941897609, 0.0412696810847696, 0.0319927106467962, 0.108815314898247, 
0.0867189417281315, 0.0472560235806836, 0.0300347697835532, 0.124252895869383, 
0.111313062025094, 0.0675699701502124, 0.0373036736422192, 0.15590486866906, 
0.153015363026187, 0.0974421665330981, 0.0563274705516883, 0.228204960789236, 
0.2051853623161, 0.169506663277659, 0.0887344403643672, 0.334252927026757, 
0.295865002768107, 0.24644310712191, 0.151125704749819), scale = c(3.33440594470802, 
3.56052439703143, 0.247924244305486, 0.150095610340916, 4.05027295845979, 
3.78800058764034, 0.35102275649495, 0.486681740716857, 3.62361340129802, 
3.48003996256557, 0.661105185961319, 0.247063080430436, 2.80272048406229, 
2.6202749242671, 0.764421647692264, 0.268131101357167, 2.29273349177322, 
2.11664063949949, 0.749125964202465, 0.417294230771058, 2.06795900869592, 
1.80352742984593, 0.78399749774511, 0.451625737966194, 1.82650979495336, 
1.62873839529095, 0.806784448546511, 0.478644431229979, 1.61927028533258, 
1.46189762327418, 0.810126410765397, 0.500531738601558, 1.53243134315106, 
1.42674067036488, 0.869134221688119, 0.569591312385404, 1.44388696689722, 
1.41174931932889, 0.938154963786188, 0.631891327069348, 1.47052388447475, 
1.4543783566239, 1.01084363142901, 0.644940081254986, 1.29365253055963, 
1.4079170674224, 1.03720449704243, 0.643101192871478, 1.08541958643012, 
1.29524623033074, 1.02176143186057, 0.673175882772333))
pop1 <- data.frame(cohort=2012:1900,
                   pop=c(rep(17239,9), 16854, 16085, 15504, 15604, 16381, 16705, 
                       16762, 16853, 15487, 14623, 14066, 13568, 13361, 13161, 13234, 
                       13088, 12472, 12142, 12062, 12078, 11426, 12027, 11963, 12435, 
                       12955, 13013, 13125, 13065, 12249, 11103, 9637, 9009, 8828, 
                       8350, 7677, 7444, 7175, 6582, 6573, 6691, 6651, 6641, 6268, 
                       6691, 6511, 6857, 7304, 7308, 7859, 7277, 8323, 8561, 7173, 
                       6942, 7128, 6819, 5037, 6798, rep(6567,46)))
## these enum strings should be moved to C++
screenT <- c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70", "fourYearlyScreen50to70", "screen50",
             "screen60", "screen70", "screenUptake", "stockholm3_goteborg",
             "stockholm3_risk_stratified", "goteborg", "risk_stratified", "mixed_screening","regular_screen","single_screen")
stateT <- c("Healthy","Localised","Metastatic")
gradeT <- c("Gleason_le_6","Gleason_7","Gleason_ge_8")
eventT <- c("toLocalised","toMetastatic","toClinicalDiagnosis",
            "toCancerDeath","toOtherDeath","toScreen","toBiopsyFollowUpScreen",
            "toScreenInitiatedBiopsy","toClinicalDiagnosticBiopsy","toScreenDiagnosis",
            "toOrganised","toTreatment","toCM","toRP","toRT","toADT","toUtilityChange","toUtility",
            "toSTHLM3", "toOpportunistic")
diagnosisT <- c("NotDiagnosed","ClinicalDiagnosis","ScreenDiagnosis")
treatmentT <- c("CM","RP","RT")
psaT <- c("PSA<3","PSA>=3") # not sure where to put this...

callFhcrc <- function(n=10,screen=screenT,nLifeHistories=10,screeningCompliance=0.75,
                      seed=12345, studyParticipation=50/260, psaThreshold=3.0, panel=FALSE,
                      includePSArecords=FALSE, flatPop = FALSE, pop = pop1, tables = IHE, debug=FALSE,
                      discountRate = 0.03, parms = NULL,
                      mc.cores=1) {
  ## save the random number state for resetting later
  state <- RNGstate(); on.exit(state$reset())
  ## yes, we use the user-defined RNG
  RNGkind("user")
  set.user.Random.seed(seed)
  ## birth cohorts that should give approximately the number of men alive in Stockholm in 2012
  ## check the input arguments
  screen <- match.arg(screen)
  stopifnot(is.na(n) || is.integer(as.integer(n)))
  stopifnot(is.integer(as.integer(nLifeHistories)))
  stopifnot(is.double(as.double(screeningCompliance)))
  screenIndex <- which(screen == screenT) - 1
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
  ## Minor changes to fhcrcData
  if (!is.null(tables)) 
      for (name  in names(tables))
          fhcrcData[[name]] <- tables[[name]]
  fhcrcData$rescreening <- rescreening
  fhcrcData$rescreening$total <- fhcrcData$rescreening$total_cat
  fhcrcData$prtx$Age <- as.double(fhcrcData$prtx$Age)
  fhcrcData$prtx$DxY <- as.double(fhcrcData$prtx$DxY)
  fhcrcData$prtx$G <- fhcrcData$prtx$G - 1L
  fhcrcData$pradt$Grade <- fhcrcData$pradt$Grade - 1L
  fhcrcData$pradt$Age <- as.double(fhcrcData$pradt$Age)
  fhcrcData$pradt$DxY <- as.double(fhcrcData$pradt$DxY)
  ## fhcrcData$biopsyComplianceTable <-
  ##     data.frame(expand.grid(psa=c(4,7,10),age=seq(55,75,by=5)),
  ##                compliance=unlist(fhcrcData$biopsy_frequency[,-(1:2),]))
  fhcrcData$biopsyOpportunisticComplianceTable <- swedenOpportunisticBiopsyCompliance 
  fhcrcData$biopsyFormalComplianceTable <- swedenFormalBiopsyCompliance
  fhcrcData$survival_local <-
      with(fhcrcData$survival_local,
           data.frame(Age=as.double(AgeLow),Grade=Grade,Time=as.double(Time),
                      Survival=Survival))
  fhcrcData$survival_dist <-
      with(fhcrcData$survival_dist,
           data.frame(Grade=Grade,Time=as.double(Time),
                      Survival=Survival))
  updateParameters <- c(parms,
                        list(nLifeHistories=as.integer(nLifeHistories),
                             screeningCompliance=screeningCompliance,
                             studyParticipation=studyParticipation,
                             psaThreshold=psaThreshold,
                             screen=as.integer(screenIndex),
                             discountRate.costs=discountRate,
                             discountRate.effectiveness=discountRate))
  parameter <- FhcrcParameters
  for (name in names(updateParameters))
      parameter[[name]] <- updateParameters[[name]]
  pind <- sapply(parameter,class)=="numeric" & sapply(parameter,length)==1
  bInd <- sapply(parameter,class)=="logical" & sapply(parameter,length)==1
  if (parameter$stockholmTreatment)
      fhcrcData$prtx <- stockholmTreatment
  ## check some parameters for sanity
  if (panel && parameter["rTPF"]>1) stop("Panel: rTPF>1 (not currently implemented)")
  if (panel && parameter["rFPF"]>1) stop("Panel: rFPF>1 (not currently implemented)")
  ## now run the chunks separately
  print(system.time(out <- parallel::mclapply(1:mc.cores,
                function(i) {
                  chunk <- chunks[[i]]
                  set.user.Random.seed(initialSeeds[[i]])
                  .Call("callFhcrc",
                        parms=list(n=as.integer(length(chunk)),
                            firstId=ns[i],
                            panel=panel, # bool
                            debug=debug, # bool
                            cohort=as.double(chunk),
                            parameter=unlist(parameter[pind]),
                            bparameter=unlist(parameter[bInd]),
                            otherParameters=parameter[!pind & !bInd],
                            tables=fhcrcData,
                            includePSArecords=includePSArecords),
                        PACKAGE="microsimulation")
                }, mc.cores = mc.cores)))
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
  ## map2df <- function(obj) "names<-"(data.frame(obj[-1]),obj[[1]]) 
  map2df <- function(obj) as.data.frame(do.call("cbind",obj))
  lifeHistories <- do.call("rbind",lapply(out,function(obj) map2df(obj$lifeHistories)))
  psarecord <- do.call("rbind",lapply(out,function(obj) data.frame(obj$psarecord)))
  falsePositives <- do.call("rbind",lapply(out,function(obj) data.frame(obj$falsePositives)))
  parameters <- map2df(out[[1]]$parameters)
  ## Identifying elements without name which also need to be rbind:ed
  costs <- do.call("rbind",lapply(out,function(obj) data.frame(obj$costs)))
  names(costs) <- c("item","cohort","age","costs")
  names(lifeHistories) <- c("id","state","ext_grade","dx","event","begin","end","year","psa")
  enum(lifeHistories$state) <- stateT
  enum(lifeHistories$dx) <- diagnosisT
  enum(lifeHistories$event) <- eventT
  enum <- list(stateT = stateT, eventT = eventT, screenT = screenT, diagnosisT = diagnosisT,
               psaT = psaT)
  out <- list(n=n,screen=screen,enum=enum,lifeHistories=lifeHistories,parameters=parameters,
              ## prev=summary$prev, pt=summary$pt, events=summary$events)
              summary=summary,costs=costs, psarecord=psarecord, cohort=data.frame(table(cohort)),
              discountRate = discountRate, falsePositives=falsePositives)
  class(out) <- "fhcrc"
  out
}

## R --slave -e "options(width=200); require(microsimulation); callFhcrc(100,nLifeHistories=1e5,screen=\"screen50\")[[\"parameters\"]]"

summary.fhcrc <- function(obj) {
    newobj <- obj[c("n","screen","discountRate")]
    with(obj,
         structure(.Data=c(newobj,
                       with(obj, list(
                           LE=sum(summary$pt$pt)/n,
                           QALE=sum(summary$ut$ut)/n,
                           costs=sum(costs$costs)/n))),
                   class="summary.fhcrc"))
}
print.summary.fhcrc <- function(obj) 
    cat(sprintf("Screening scenario: \t%s
Life expectancy: \t%f
Discounted QALE: \t%f
Discounted costs: \t%f
Discounted rate: \t%f
",obj$screen,obj$LE,obj$QALE,obj$costs,obj$discountRate))

ICER <- function(object1, object2, ...)
    UseMethod("ICER")

ICER.fhcrc <- function(obj1,obj2,...) {
    stopifnot(obj1$discountRate == obj2$discountRate)
    summary1 <- summary(obj1,...)
    summary2 <- summary(obj2,...)
    out <- list(ICER.QALE=(summary1$costs-summary2$costs)/(summary1$QALE-summary2$QALE),
                delta.QALE=summary1$QALE-summary2$QALE,
                delta.costs=summary1$costs-summary2$costs)
    if (obj1$discountRate == 0)
        out <- c(out,
                 list(ICER.LE = (summary1$costs-summary2$costs)/(summary1$LE-summary2$LE),
                      delta.LE = summary1$LE-summary2$LE))
    out
}

print.fhcrc <- function(obj,...)
    cat(sprintf("FHCRC prostate cancer model with %i individual(s) under scenario '%s'.\n",
                obj$n, obj$screen),
        ...)

plot.fhcrc <- function(obj,type=c("incidence","cancerdeath"),plot.type="l",xlim=c(40,100), add=FALSE, ...) {
    type <- match.arg(type)
    event_types <- switch(type,
                          incidence=c("toClinicalDiagnosis","toScreenDiagnosis"),
                          cancerdeath="toCancerDeath")
    if (require(dplyr)) {
        pt <- obj$summary$pt %>%
            group_by(age) %>%
                summarise(pt=sum(pt))
        events <- obj$summary$events %>%
            filter(event %in% event_types) %>%
                group_by(age) %>%
                    summarise(n=sum(n))
        out <- left_join(pt,events,by="age") %>% mutate(rate = ifelse(is.na(n), 0, n/pt))
        if (!add) plot(rate~age, data=out, type=plot.type, xlim=xlim, ...) else lines(rate~age, data=out,  ...)
    } else error("dplyr is not available for plotting")
}

lines.fhcrc <- function(obj,...) {
    plot(obj, ..., add=TRUE)
}



## utility - not exported
assignList <- function(lst,...)
  for(i in 1:length(lst))
    assign(names(lst)[i], lst[[i]], ...)
## assignList(formals(callFhcrc),pos=1)

NN.fhcrc <- function(obj, ref.obj, startAge = 50, stopAge = Inf) {
    if (require(dplyr)) {
        pNNS <- function(thisScenario) { 
            as.numeric((thisScenario$summary$events %>%
                        filter(event=="toCancerDeath" & age>=startAge & age<stopAge) %>%
                        summarise(sumEvents=sum(n))) / # divided by
                       (thisScenario$summary$prev %>%
                        filter(age==round(startAge)) %>%
                        summarise(sumPop=sum(count))))
        }
        pNND <- function(thisScenario) {
            as.numeric(thisScenario$summary$events %>%
                       filter(event=="toCancerDeath" & age>=startAge & age<stopAge) %>%
                       summarise(sumEvents=sum(n)) /
                       thisScenario$summary$events %>%
                       filter(is.element(event,c("toScreenDiagnosis","toClinicalDiagnosis")) & age>=startAge & age<stopAge) %>%
                       summarise(sumEvents=sum(n)))
        }
        NNS <- 1 / (pNNS(ref.obj) - pNNS(obj)) #number needed to screen to prevent 1 PCa death
        NND <- 1 / (pNND(ref.obj) - pNND(obj)) #number needed to detect to prevent 1 PCa death
        ## Include additional number needed to treat (NNT) [Gulati 2011] to show overdiagnosis?
        return(list(NNS=NNS,NND=NND))
    } else error("NN.fhcrc: require dplyr to calculate NNS and NND")
}

ggplot.fhcrc <- function(obj,type=c("psa","biopsies","incidence","metastatic","cancerdeath","alldeath"),ages=c(50,85), ...) {
    type <- match.arg(type)
    event_types <- switch(type,
                          psa="toScreen",
                          biopsies=c("toClinicalDiagnosticBiopsy","toScreenInitiatedBiopsy"),
                          incidence=c("toClinicalDiagnosis","toScreenDiagnosis"),
                          metastatic="toMetastatic",
                          cancerdeath="toCancerDeath",
                          alldeath=c("toCancerDeath","toOtherDeath"))
    if(class(obj)!="list"){obj <- list(obj)}
    if (require(ggplot2) & require(dplyr)) {
        pt <- do.call("rbind",lapply(obj,function(obj) cbind(obj$summary$pt,pattern=obj$screen))) %>%
            group_by(pattern,age) %>%
                summarise(pt=sum(pt))
        events <-  do.call("rbind",lapply(obj,function(obj) cbind(obj$summary$events,pattern=obj$screen))) %>%
            filter(is.element(event, event_types)) %>%
                group_by(pattern,age) %>%
                    summarise(n=sum(n))
        out <- left_join(pt,events,by=c("pattern","age")) %>%
            mutate(rate = 1000*ifelse(is.na(n), 0, n/pt)) %>%
                filter(age >= min(ages),
                       age <= max(ages))
        ggplot(out, aes(age, rate, group=pattern, colour=pattern)) + ...
    } else error("ggplot.fhcrc: require both ggplot2 and dplyr")
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
