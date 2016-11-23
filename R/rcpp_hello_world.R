
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
      if (length(ichull)<2) return(ichull)
      xi <- x[ichull]
      yi <- y[ichull]          # subset to convex hull
      imin <- which(xi==min(xi))
      include <- sapply(1:length(ichull),
                        function(i)       # establish the frontier
                        i==imin || (i>imin && yi[i-1]<yi[i] && xi[i-1]<xi[i]))
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
    revised_natural_history=TRUE,
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
    mubeta2=c(0.0397,0.1678, 0.0), # base::grade
    sebeta2=c(0.0913,0.3968, 0.0), # base::grade
    rev_mubeta2=c(0.051, 0.129, 0.1678), # ext::grade
    rev_sebeta2=c(0.064, 0.087, 0.3968), # ext::grade
    alpha7=0.01241490,
    beta7=0.01543417,
    alpha8=-4.60517081,
    beta8=0.08770374,
    gamma_m_diff=0.0001,
    RR_T3plus=2.0,
    ## mubeta2.scale=1.0, # cf. 2.1
    ## beta.rho=0.62,
    c_txlt_interaction = 1.0,
    c_baseline_specific = 1.0,
    c_benefit_value0 = 10, # (value -> reduction): 0.04 -> 10%; 0.18 -> 20%; 10 -> 28%
    sxbenefit = 1.0,
    c_benefit_type = 0, # 0=stage-shift (=> c_benefit_value0=10), 1=lead-time based (=> c_benefit_value1=0.1)
    c_benefit_value1 = 0.1,
    screeningCompliance = 0.75,
    rescreeningCompliance = 0.95,
    biopsyCompliance = 0.858,
    biopsySensitivity = 0.8,
    studyParticipation = 50.0/260.0,
    nLifeHistories = 10L, screen = 0L, ## integers
    psaThreshold = 3.0,
    psaThresholdBiopsyFollowUp = 4.0,
    ## BPThreshold=3.47,
    ## BPThresholdBiopsyFollowUp=3.47,
    ## BPThreshold=4.69,
    ## BPThresholdBiopsyFollowUp=4.69,
    biomarker_model = 0, # biomarker_model = 0 random, biomarker_model = 1 psa/risk based correction of FP
    PSA_FP_threshold_nCa=4.15, # reduce FP in no cancers with PSA threshold
    PSA_FP_threshold_GG6=3.41, # reduce FP in GG 6 with PSA threshold
    BPThreshold=4.2,
    BPThresholdBiopsyFollowUp=4.2,
    ## Natural history calibration
    gleason_le_6_hr = 1,
    gleason_7_hr = 1,
    gleason_ge_8_hr = 1,
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
    hr_locoregional=transform(expand.grid(age=c(50,60,70),ext_grade=0:2,psa10=0:1),
                              hr=c(0.4852788, 0.7309101, 1.5703646,
                                   3.1780491, 2.2476271, 3.1391895,
                                   1.2743633, 0.9972281, 1.2596770,
                                   1.0367029, 0.9528743, 1.3431982,
                                   7.3491808, 3.4922388, 2.8214167,
                                   0.7880672, 0.7624732, 0.7827559)),
    hr_metastatic=data.frame(age=c(50, 60, 70),
                             hr=c(0.8325735, 0.9403021, 0.7998358)),
    cost_parameters = c("Invitation" = 50,
                        "Formal PSA" = 130,
                        "Formal panel" = 730,
                        "Opportunistic PSA" = 1910,
                        "Opportunistic panel" = 2510, #N.B. This one is new and should be used
                        "Biopsy" = 12348,
                        "Prostatectomy" = 117171,
                        "Radiation therapy" = 117171,
                        "Active surveillance" = 141358,
                        "Cancer death" = 585054),
    ## IHE doesn't use the postrecovery period (as reported in the Heijnsdijk 2012 reference), should we?
    production = data.frame(ages = c(0, 55, 65, 75),
                            values=c(467433.137375286, 369392.309986899, 45759.6141748681, 0.0)),
    lost_production_proportions= c("Formal PSA"=0.0011,
                                   "Formal panel"=0.0011,
                                   "Opportunistic PSA"=0.0025,
                                   "Opportunistic panel"=0.0025,
                                   "Biopsy"=0.0044,
                                   "Prostatectomy"=0.1083,
                                   "Radiation therapy"=0.1250,
                                   "Active surveillance"=0.0833,
                                   "Metastatic cancer"=0.7602),
    utility_estimates = 1 - c("Invitation" = 1,
                              "Formal PSA" = 0.99,
                              "Formal panel" = 0.99,
                              "Opportunistic PSA" = 0.99,
                              "Opportunistic panel" = 0.99,
                              "Biopsy" = 0.90,
                              "Prostatectomy part 1" = 0.67,
                              "Prostatectomy part 2" = 0.77,
                              "Radiation therapy part 1" = 0.73,
                              "Radiation therapy part 2" = 0.78,
                              "Active surveillance" = 0.97,
                              "Palliative therapy" = 0.60,
                              "Terminal illness" = 0.40,
                              "Metastatic cancer" = 0.85,
                              "Death" = 0.00),
    ## Utility duration is given in years.
    utility_duration = c("Invitation" = 0.0,
                         "Formal PSA" = 1/52,
                         "Formal panel" = 1/52,
                         "Opportunistic PSA" = 1/52,
                         "Opportunistic panel" = 1/52,
                         "Biopsy" = 3/52,
                         "Prostatectomy part 1" = 2/12,
                         "Prostatectomy part 2" = 10/12,
                         "Radiation therapy part 1" = 2/12,
                         "Radiation therapy part 2" = 10/12,
                         "Active surveillance" = 7,
                         "Palliative therapy" = 30/12,
                         "Terminal illness" = 6/12)
)
IHE <- list(prtx=data.frame(Age=50.0,DxY=1973.0,G=1:2,CM=0.6,RP=0.26,RT=0.14)) ## assumed constant across ages and periods
ParameterNV <- FhcrcParameters[sapply(FhcrcParameters,class)=="numeric" & sapply(FhcrcParameters,length)==1]
## ParameterIV <- FhcrcParameters[sapply(FhcrcParameters,class)=="integer" & sapply(FhcrcParameters,length)==1]
swedenOpportunisticBiopsyCompliance <- data.frame(
    psa = c(3, 5, 10, 3, 5, 10, 3, 5, 10, 3, 5, 10, 3, 5, 10),
    age = c(40, 40, 40, 50, 50, 50, 60, 60, 60, 70, 70, 70, 80, 80, 80),
    compliance = c(0.3764045, 0.5680751, 0.7727273, 0.3110770, 0.5726548, 0.7537372, 0.2385155, 0.4814588, 0.6929770, 0.1754264, 0.3685056, 0.5602030, 0.1629213, 0.2697368, 0.5010052))
swedenFormalBiopsyCompliance <- cbind(expand.grid(psa=c(3,5,10),age=seq(40,80,10)),
                                compliance=0.858)
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
ext_stateT <- c("Healthy","T1_T2","T3plus","Metastatic")
gradeT <- c("Gleason_le_6","Gleason_7","Gleason_ge_8","Healthy")
eventT <- c("toLocalised","toMetastatic","toClinicalDiagnosis",
            "toCancerDeath","toOtherDeath","toScreen","toBiopsyFollowUpScreen",
            "toScreenInitiatedBiopsy","toClinicalDiagnosticBiopsy","toScreenDiagnosis",
            "toOrganised","toTreatment","toCM","toRP","toRT","toADT","toUtilityChange","toUtility",
            "toSTHLM3", "toOpportunistic")
diagnosisT <- c("NotDiagnosed","ClinicalDiagnosis","ScreenDiagnosis")
treatmentT <- c("no_treatment","CM","RP","RT")
psaT <- c("PSA<3","PSA>=3") # not sure where to put this...

callFhcrc <- function(n=10,screen=screenT,nLifeHistories=10,
                      seed=12345,
                      panel=FALSE,
                      includePSArecords=FALSE, includeDiagnoses=FALSE,
                      flatPop = FALSE, pop = pop1, tables = IHE, debug=FALSE,
                      parms = NULL, mc.cores=1, ...) {
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
                             screen=as.integer(screenIndex)))
  parameter <- FhcrcParameters
    for (name in names(updateParameters)) {
        if(!(name %in% names(parameter))
           warning("Name in parms argument not in FhcrcParameters: ",name,".")
      parameter[[name]] <- updateParameters[[name]]
      }
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
                            includePSArecords=includePSArecords,
                            includeDiagnoses=includeDiagnoses),
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
  diagnoses <- do.call("rbind",lapply(out,function(obj) data.frame(obj$diagnoses)))
  falsePositives <- do.call("rbind",lapply(out,function(obj) data.frame(obj$falsePositives)))
  parameters <- map2df(out[[1]]$parameters)
  ## Identifying elements without name which also need to be rbind:ed
  societal.costs <- do.call("rbind",lapply(out,function(obj) data.frame(obj$costs))) #split in sociatal and healthcare perspective
  ## names(costs) <- c("type","item","cohort","age","costs")
  names(societal.costs) <- c("type","item","age","costs")
  societal.costs$type <- ifelse(societal.costs$type, "Productivity loss", "Health sector cost") # societal perspective
  healthsector.costs <- societal.costs[societal.costs["type"] == "Health sector cost", c("item", "age", "costs")] # healthcare perspective
  names(lifeHistories) <- c("id", "state", "ext_grade", "dx", "event", "begin", "end", "year", "psa", "utility")
  enum(lifeHistories$state) <- stateT
  enum(lifeHistories$dx) <- diagnosisT
  enum(lifeHistories$event) <- eventT
  enum(diagnoses$state) <- stateT
  enum(diagnoses$ext_state) <- ext_stateT
  enum(diagnoses$ext_grade) <- gradeT
  enum(diagnoses$dx) <- diagnosisT
  enum(diagnoses$tx) <- treatmentT
  enum <- list(stateT = stateT, eventT = eventT, screenT = screenT, diagnosisT = diagnosisT,
               psaT = psaT, ext_stateT = ext_stateT)
  out <- list(n=n,screen=screen,enum=enum,lifeHistories=lifeHistories,
              parameters=parameters, summary=summary,
              healthsector.costs=healthsector.costs, societal.costs=societal.costs,
              psarecord=psarecord, diagnoses=diagnoses,
              cohort=data.frame(table(cohort)),simulation.parameters=parameter,
              falsePositives=falsePositives)
  class(out) <- "fhcrc"
  out
}

## R --slave -e "options(width=200); require(microsimulation); callFhcrc(100,nLifeHistories=1e5,screen=\"screen50\")[[\"parameters\"]]"

summary.fhcrc <- function(object, ...) {
    newobj <- object[c("n","screen")]
    with(object,
         structure(.Data=c(newobj,
                       with(object, list(
                           discountRate.costs=simulation.parameters$discountRate.costs,
                           discountRate.effectiveness=simulation.parameters$discountRate.effectiveness,
                           LE=sum(summary$pt$pt)/n,
                           QALE=sum(summary$ut$ut)/n,
                           healthsector.costs=sum(healthsector.costs$costs)/n,
                           societal.costs=sum(societal.costs$costs)/n))),
                   class="summary.fhcrc"))
}

print.summary.fhcrc <- function(x, ...) {
    obj <- x
    cat(sprintf(
"Screening scenario:            %s
Life expectancy:                %f
Discounted QALE:                %f
Discounted health sector costs: %f
Discounted societal costs:      %f
Discounted rate (effect.):      %f
Discounted rate (costs):        %f
", obj$screen, obj$LE, obj$QALE,
obj$healthsector.costs, obj$societal.costs,
obj$discountRate.effectiveness,
obj$discountRate.costs))
}

ICER <- function(object1, object2, ...)
    UseMethod("ICER")

ICER.fhcrc <- function(object1, object2,
                       perspective = c("societal.costs",
                                       "healthsector.costs"), ...) {
    perspective <- match.arg(perspective)
    p1 <- object1$simulation.parameters
    p2 <- object2$simulation.parameters
    stopifnot(p1$discountRate.costs == p2$discountRate.costs)
    stopifnot(p1$discountRate.effectiveness == p2$discountRate.effectiveness)
    summary1 <- summary(object1, ...)
    summary2 <- summary(object2, ...)
    out <- list(ICER.QALE=(summary1[[perspective]] - summary2[[perspective]]) /
                    (summary1$QALE - summary2$QALE),
                delta.QALE=summary1$QALE - summary2$QALE,
                delta.costs=summary1[[perspective]] - summary2[[perspective]])
    if (p1$discountRate.costs == 0 && p2$discountRate.costs == 0 &&
        p1$discountRate.effectiveness == 0 && p2$discountRate.effectiveness == 0) {
        out <- c(out,
                 list(ICER.LE = (summary1[[perspective]] - summary2[[perspective]]) /
                          (summary1$LE - summary2$LE),
                      delta.LE = summary1$LE - summary2$LE))
    }
    out
}

print.fhcrc <- function(x, ...)
    cat(sprintf("FHCRC prostate cancer model with %i individual(s) under scenario '%s'.\n",
                x$n, x$screen),
        ...)

## fast operations by group using base-R
## http://stackoverflow.com/questions/3685492/r-speeding-up-group-by-operations
grp_apply = function(XS, INDEX, FUN, ..., simplify=T) {
  FUN = match.fun(FUN)
  if (!is.list(XS))
    XS = list(XS)
  as.data.frame(as.table(tapply(1:length(XS[[1L]]), INDEX, function(s, ...)
    do.call(FUN, c(lapply(XS, `[`, s), list(...))), ..., simplify=simplify)))
}

## TODO: include prevalences and relative rate-ratios in switch. Also
## allow for ceiling on groups to allow for other than yearly rates
## for the time
predict.fhcrc <- function(object, scenarios=NULL,
                          type = "incidence.rate", group = "age", ...) {
    if(!inherits(object,"fhcrc")) stop("Expecting object to be an fhcrc object")
    if(!(is.null(scenarios) || all(sapply(scenarios,inherits,"fhcrc")) || inherits(object,"fhcrc")))
        stop("Expecting scenarios is NULL, a fhcrc object or a list of fhcrc objects")

    ## Stripping of potential rate ratio option before matching
    abbr_type <- match.arg(sub(".?rr$|.?rate.?ratio$",
                               "", type, ignore.case = TRUE),
                           c("incidence.rate", "testing.rate",
                             "biopsy.rate", "metastasis.rate",
                             "pc.mortality.rate",
                             "allcause.mortality.rate", "prevalence"))

    ## Allowing for several groups
    group <- match.arg(group,
                       c("state", "grade", "dx", "psa", "age", "year"),
                       several.ok = TRUE)

    event_types <- switch(abbr_type,
                          incidence.rate = c("toClinicalDiagnosis", "toScreenDiagnosis"),
                          testing.rate = "toScreen",
                          biopsy.rate = c("toClinicalDiagnosticBiopsy", "toScreenInitiatedBiopsy"),
                          metastasis.rate = "toMetastatic",
                          pc.mortality.rate = "toCancerDeath",
                          allcause.mortality.rate = c("toCancerDeath", "toOtherDeath"))

    ## Fixes colnames after group operation
    name_grp <- function(x) {names(x)[grep("^Var[0-9]+$", names(x))] <- group; x}

    ## Calculates rates of specific events by specified groups
    calc_rate <- function(object, event_types, group){
        pt <- with(object$summary$pt,
                   name_grp(grp_apply(pt,
                                      lapply(as.list(group), function(x) eval(parse(text = x))),
                                      sum)))
        ## temp fix: no events causes angst
        ## todo: if subset has no dim replace with zeros
        if(!any(object$summary$event$event %in% event_types)) {
            stop(paste("The event(s)", paste(event_types, collapse = ", "),
                       "was not found in the", object$screen, "scenario"))
        }
        events <- with(subset(object$summary$events, event %in% event_types),
                       name_grp(grp_apply(n,
                                          lapply(as.list(group), function(x) eval(parse(text = x))),
                                          sum)))
        within(merge(pt, events, by = group, all = TRUE),{
            if("age" %in% group) age <- as.numeric(levels(age))[age] #important factor conversion
            if("year" %in% group) year <- as.numeric(levels(year))[year] #important factor conversion
            rate <- ifelse(is.na(Freq.y) & !is.na(Freq.x), 0, Freq.y/Freq.x) #no events but some pt -> 0
            n <- Freq.y
            pt <- Freq.x
            rm(Freq.x,Freq.y)})
    }

    ## Calculate prevalences by specified groups
    calc_prev <- function(object, group){
        within(with(object$summary$prev,
                    name_grp(
                        grp_apply(count,
                                  lapply(as.list(group),
                                         function(x) eval(parse(text = x))), sum))),{
                                             if("age" %in% group) age <- as.numeric(levels(age))[age] #important factor conversion
                                             if("year" %in% group) year <- as.numeric(levels(year))[year] #important factor conversion
                                             prevalence <- Freq/object$n
                                             rm("Freq")
                                         })
    }

    ## Calculates the outcome in the passed function for all
    ## simulation objects in the 'scenarios' list. Then the object
    ## outcome (e.g. rates or prev) are for the scenarios are added as
    ## rows and the scenario name as a column.
    predict_scenarios <- function(scenarios, calc_outcome, ...) {
        do.call(rbind, lapply(scenarios,
        {function(object, ...)
            cbind(calc_outcome(object, ...), scenario = object$screen)}, ...))
    }

    ## Input checks allow for scenarios to be a single fhcrc object or
    ## list of fhcrc objects. Now make sure scenarios is a list.
    if(inherits(scenarios, "fhcrc")) {scenarios <- list(scenarios)}

    ## Rate-ratio if type ends with rate.ratio or RR
    if(grepl(".?rate.?ratio$|.?rr$", type, ignore.case = TRUE)){
        scenario_rates <- predict_scenarios(unique(scenarios), calc_rate, event_types, group)
        reference_rate <- predict_scenarios(list(object), calc_rate, event_types, group)
        within(merge(scenario_rates, reference_rate, by = group),{
            scenario <- scenario.x
            rate.ratio <- rate.x/rate.y
            rate.ratio[!is.finite(rate.ratio)] <- NaN
            rm(list=ls(pattern=".x$|.y$"))})

        ## Prevalence if type ends with rate.ratio or RR
    } else if(grepl(".?prev$|.?prevalence$", type, ignore.case = TRUE)){
        predict_scenarios(unique(c(list(object),scenarios)), calc_prev, group)

        ## Defauls to plain rates. If reference object exist add it to
        ## scenario list and remove duplicates.
    }else{
        predict_scenarios(unique(c(list(object),scenarios)), calc_rate, event_types, group)
    }
}

plot.fhcrc <- function(x, type=c("incidence.rate", "testing.rate",
                                 "biopsy.rate", "metastasis.rate",
                                 "pc.mortality.rate",
                                 "allcause.mortality.rate"),
                       plot.type="l", add=FALSE, xlab="Age (years)",
                       ylab=NULL, ...) {
    type <- match.arg(type)
    if (is.null(ylab)) {ylab <- switch(type,
                                       incidence.rate="Prostate cancer incidence rates per 100,000",
                                       testing.rate="PSA rates per 1000",
                                       biopsy.rate="Biopsies per 1000",
                                       metastasis.rate="Metastatic onset per 100,000",
                                       pc.mortality.rate="Cancer mortality rates per 100,000",
                                       allcause.mortality.rate="All cause mortality rates per 100,000")}
    rates <- predict(object = x, type = type)
    rates$rate = rates$rate*switch(type, testing.rate=1000, biopsy.rate=1000, incidence.rate=1e5, metastasis.rate=1e5,pc.mortality.rate=1e5,allcause.mortality.rate=1e5)
    if (!add) plot(rate~age, data=rates, type=plot.type, xlab=xlab, ylab=ylab, ...) else lines(rate~age, data=rates,  ...)
}
lines.fhcrc <- function(x,...) {
    plot(x, ..., add=TRUE)
}

## utility - not exported
assignList <- function(lst,...)
  for(i in 1:length(lst))
    assign(names(lst)[i], lst[[i]], ...)
## assignList(formals(callFhcrc),pos=1)

NN.fhcrc <- function(obj, ref.obj, startAge = 50, stopAge = Inf) {
    pNNS <- function(thisScenario) {
        with(subset(thisScenario$summary$events,
                    event=="toCancerDeath" & age>=startAge & age<stopAge),
             sum(n)) / # divided by
            with(subset(thisScenario$summary$prev,
                        age==round(startAge)),
                 sum(count))
    }
    pNND <- function(thisScenario) {
        with(subset(thisScenario$summary$events,
                    event=="toCancerDeath" & age>=startAge & age<stopAge),
             sum(n)) / # divided by
            with(subset(thisScenario$summary$events,
                        event %in% c("toScreenDiagnosis","toClinicalDiagnosis") & age>=startAge & age<stopAge),
                 sum(n))
    }
    NNS <- 1 / (pNNS(ref.obj) - pNNS(obj)) #number needed to screen to prevent 1 PCa death
    NND <- 1 / (pNND(ref.obj) - pNND(obj)) #number needed to detect to prevent 1 PCa death
    ## Include additional number needed to treat (NNT) [Gulati 2011] to show overdiagnosis?
    return(list(NNS=NNS,NND=NND))
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
                    help = function() {
                        'Reference class implementation of an event queue. Fields for the event times and the events list.'
                    },
                    insert = function(time,event) {
                        'Method to insert the event at the given time'
                        insert.ord <- findInterval(time,times)
                        times <<- append(times,time,insert.ord)
                        events <<- append(events,event,insert.ord)
                    },
                    pop = function() {
                        'Method to remove the head of the event queue and return its value'
                        head <- structure(events[[1]], time=times[1])
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
                    remove = function(predicate, ...) {
                        'Method to remove events that satisfy some predicate'
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
                    help = function() {
                        'Reference class implementation of an event-oriented discrete event simulation. Fields for the event times and the events list.'
                    },
                    scheduleAt = function(time, event) {
                        'Method that adds attributes for the event time and the sendingTime, and then insert the event into the event queue'
                        attr(event,"time") <- time
                        attr(event,"sendingTime") <- currentTime
                        insert(time, event)
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
