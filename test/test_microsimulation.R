## try(detach("package:microsimulation", unload=TRUE))
## require(microsimulation)
## microsimulation:::.testPackage()

## unit tests
require(microsimulation)
temp=callCalibrationPerson(10)
stopifnot(temp$StateOccupancy[1:2] == c(422,354))
temp2=callFhcrc(1000)
stopifnot(abs(with(temp2,sum(summary$pt$pt)/n)-79.92847)<1e-3)
temp3 <- callIllnessDeath(10)
stopifnot(abs(with(temp3,sum(pt$pt)/10)-64.96217)<1e-3)

temp=callCalibrationPerson(10)
stopifnot(temp$StateOccupancy[1:2] == c(422,354))
temp3 <- callIllnessDeath(10)
stopifnot(abs(with(temp3,sum(pt$pt)/10)-64.96217)<1e-3)


## all(c(callFhcrc(5)$lifeHistories == callFhcrc(5)$lifeHistories,
##       callFhcrc(5)$parameters == callFhcrc(5)$parameters,
##       callFhcrc(5)$summary$pt == callFhcrc(5)$summary$pt,
##       callFhcrc(5)$summary$prev == callFhcrc(5)$summary$prev))
## all(callFhcrc(10)$lifeHistories == callFhcrc(10)$lifeHistories) # fails
## all(callFhcrc(1e4)$parameters == callFhcrc(1e4)$parameters) # okay

## extract PSA values
refresh
require(microsimulation)
set.seed(12345)
temp <- callFhcrc(350000,screen="stockholm3_risk_stratified",includePSArecords=TRUE)$psarecord
temp2 <- subset(temp,organised & age>=50 & age<70 & !dx) 
## first organised screen
i <- tapply(1:nrow(temp2),temp2$id,min)
temp2 <- temp2[i,]
xtabs(~state+ext_grade+I(psa>=3),temp2)
temp2 <- transform(temp2, advanced=(state>0 & ext_grade==2), cancer=(state>0))
tau2 = 0.0829 # variance

set.seed(12345+1)
temp2 <- transform(temp2,
                   psa2=Z+0.3+rnorm(nrow(temp2),0,sqrt(tau2)),
                   BBP=Z+rnorm(nrow(temp2),0,0.5*sqrt(tau2)))
## STHLM3 simulation report
with(list(threshold=3),with(temp2,
     cat(sprintf("
Advanced:\t\t\t%i
Cancers:\t\t\t%i
PSA>=3 & advanced:\t\t%i
PSA>=3 & cancer:\t\t%i
PSA>=1 & BBP>=%3.1f & adv:\t%i
PSA>=1 & BBP>=%3.1f & can:\t%i
PSA>=1 & BBP>=%3.1f & PSA<3:\t%i
BBP<%3.1f & PSA>=3:\t\t%i
(PSA>=1 & BBP>=%3.1f) | PSA>=3:\t%i
Prop reduction in biospies:\t%5.3f\n",
                 sum(advanced),
                 sum(cancer),
                 sum(psa2>=3 & advanced),
                 sum(psa2>=3 & cancer),
                 threshold,sum(psa2>=1 & BBP>=threshold & advanced),
                 threshold,sum(psa2>=1 & BBP>=threshold & cancer),
                 threshold,sum(psa2>=1 & BBP>=threshold & psa2<3),
                 threshold,sum(psa2>=1 & BBP<threshold & psa2>=3),
                 threshold,sum(psa2>=1 & BBP>=threshold & psa2>=3),
                 (sum(psa2>=1 & BBP<threshold & psa2>=3) - sum(psa2>=1 & BBP>=threshold & psa2<3))/
                       sum(psa2>=1 & BBP>=threshold & psa2>=3)
))))
## The baseline FHCRC model assumes that PSA is an unbiased measure of the underlying diease process. The results here suggest that imprecision in the measure is less important than bias - and that PSA would need to be relatively biased to get the predicted change in biopsies from STHLM3.

## The main challenge now is that the FHCRC model was based on PCPT trial data which will not be available - nor, probably, will the bias be estimable from observed data.

## testing the user-defined random number generator
init.seed <- as.integer(c(407,rep(12345,6)))
RNGkind("user")
set.user.Random.seed(init.seed)
testA <- runif(2)
next.user.Random.substream()
testB <- runif(2)
set.user.Random.seed(parallel::nextRNGStream(init.seed))
newSeed <- user.Random.seed()
testC <- runif(2)
set.user.Random.seed(parallel::nextRNGStream(newSeed))
testD <- runif(2)
##
RNGkind("L'Ecuyer-CMRG")
init.seed <- as.integer(c(407,rep(12345,6)))
.Random.seed <- init.seed
all(testA == runif(2))
.Random.seed <- parallel::nextRNGSubStream(init.seed)
all(testB == runif(2))
newSeed <- .Random.seed <- parallel::nextRNGStream(init.seed)
all(testC == runif(2))
.Random.seed <- parallel::nextRNGStream(newSeed)
all(testD == runif(2))

## More unit tests required

system.time(callFhcrc(1e5))
system.time(callFhcrc(1e5,mc.cores=4))
system.time(callFhcrc(1e6,mc.cores=4))

## Reading in the data from FHCRC
fhcrcData <- lapply(dir("~/src/fhcrc/data")[-10],
               function(name) structure(read.table(paste("~/src/fhcrc/data/",
                                                         name,sep=""),
                                                   head=TRUE,sep=","),
                                        filename=name))
lookup <- data.frame(filename=c("all_cause_mortality.csv",
   "biopsy_frequency.csv",
   "biopsy_sensitivity_smoothed.csv",
   "seer_incidence_imputed.csv",
   "hormone_frequency.csv",
   "primary_treatment_frequency.csv",
   "dre_sensitivity.csv",
   "gleason_7_frequency.csv",
   "stage_T2a_frequency.csv",
   "prostate_cancer_survival_local-regional.csv",
   "prostate_cancer_survival_distant.csv"),
  enum=c("all_cause_mortality",
    "biopsy_frequency",
    "biopsy_sensitivity",
    "obs_incidence",
    "pradt",
    "prtx",
    "dre",
    "prob_grade7",
    "prob_earlystage",
    "survival_local",
    "survival_dist"), stringsAsFactors = FALSE)
lookup <- subset(lookup, enum!="obs_incidence")
names(fhcrcData) <- with(lookup, enum[match(lapply(fhcrcData,attr,"filename"),
                                            filename)])
save("fhcrcData",file="~/src/R/microsimulation/data/fhcrcData.rda")
## lapply(fhcrcData,head)
##
## biopsy frequency
## with(fhcrcData[[2]],data.frame(psa=rep(PSA.beg,5),
##                           age=rep(c(55,60,65,70,75),each=3),
##                           biopsy_frequency=unlist(temp[[2]][,-(1:2)])))
## fhcrcData[[2]]

## testing using parallel
require(parallel)
require(microsimulation)
n <- 1e4
system.time(test <- mclapply(1:10,
                             function(i) callFhcrc(n,screen="noScreening"),
                             mc.cores=1))
system.time(test <- mclapply(1:10,
                             function(i) callFhcrc(n,screen="noScreening"),
                             mc.cores=4))
##
test <- lapply(1:10, function(i) callFhcrc(10,screen="noScreening"))
test2 <- list(lifeHistories=do.call("rbind", lapply(test,function(obj) obj$lifeHistories)),
              enum=test[[1]]$enum,
              n=sum(sapply(test,function(obj) obj$n)),
              parameters=do.call("rbind", lapply(test,"[[", "parameters")),
              summary=list(pt=do.call("rbind", lapply(test,function(obj) obj$summary$pt)),
                events=do.call("rbind", lapply(test,function(obj) obj$summary$events)),
                prev=do.call("rbind", lapply(test,function(obj) obj$summary$prev))))

## baseline analysis
options(width=110)
require(microsimulation)
n <- 1e7
n.cores <- 4
compliance <- 0.75
participation <- 1.0
noScreening <- callFhcrc(n,screen="noScreening",mc.cores=n.cores)
## "screenUptake", "stockholm3_goteborg", "stockholm3_risk_stratified"
uptake <- callFhcrc(n,screen="screenUptake",mc.cores=n.cores,
                    studyParticipation=participation,
                    screeningCompliance=compliance)
goteborg <- callFhcrc(n,screen="stockholm3_goteborg",mc.cores=n.cores,
                      studyParticipation=participation,
                      screeningCompliance=compliance)
riskStrat <- callFhcrc(n,screen="stockholm3_risk_stratified",mc.cores=n.cores,
                       studyParticipation=participation,
                       screeningCompliance=compliance)

## Lexis diagrams
plotLexis <- function(obj) {
    stopifnot(require(Epi))
    stopifnot(require(sqldf))
    history <- obj$lifeHistories
    param <- obj$parameters
    tab <- sqldf("select t1.*, ageAtCancerDiagnosis, cohort, t0 from (select id, end as ageAtDeath, (event='toCancerDeath') as cancerDeath from history where event in ('toOtherDeath','toCancerDeath')) as t1  inner join param as p on p.id=t1.id left join (select id, end as ageAtCancerDiagnosis from history where event in ('toClinicalDiagnosis','toScreenDiagnosis')) as t2 on t1.id=t2.id")
    lexis1 <- Lexis(entry=list(coh=cohort,age=0),exit=list(coh=cohort+ageAtDeath,age=ageAtDeath),
                    data=tab)
    plot(lexis1, xlab="Calendar period", ylab="Age (year)", ylim=c(0,100), asp=1)
    with(subset(tab,!is.na(ageAtCancerDiagnosis)),
         points(cohort+ageAtCancerDiagnosis,ageAtCancerDiagnosis,pch=19,cex=0.4,col="red"))
    with(subset(tab,t0+35<ageAtDeath),
         points(cohort+t0+35,t0+35,pch=19,cex=0.2,col="blue"))
    legend("topleft",legend=c("Latent cancer onset","Cancer diagnosis"),
           pch=19,col=c("blue","red"),bty="n")
}
pdf(file="~/work/lexis-20131128.pdf",width=5,height=4)
par(mar=c(5.1, 4.1, 4.1-2, 2.1))
plotLexis(noScreening)
dev.off()

## rate calculations
pop <- data.frame(age=0:100,pop=c(12589, 14785, 15373, 14899, 14667,
14437, 14076, 13386, 13425, 12971, 12366, 11659, 11383, 10913, 11059,
11040, 11429, 12303, 13368, 13388, 13670, 13539, 13886, 13913, 14269,
14508, 15073, 15419, 15767, 15721, 16328, 16489, 17126, 16345, 15573,
15702, 16017, 16251, 17069, 16853, 16898, 16506, 15738, 15151, 15224,
15960, 16248, 16272, 16325, 14963, 14091, 13514, 13000, 12758, 12521,
12534, 12333, 11699, 11320, 11167, 11106, 10427, 10889, 10732, 11042,
11367, 11269, 11210, 10982, 10115, 9000, 7652, 6995, 6680, 6144, 5473,
5108, 4721, 4130, 3911, 3756, 3507, 3249, 2803, 2708, 2355, 2188,
2020, 1734, 1558, 1183, 1064, 847, 539, 381, 277, 185, 90, 79, 48,
61))
w <- with(subset(pop,age>=50 & age<80),data.frame(age=age,wt=pop/sum(pop)))

require(sqldf)
eventRates <- function(obj,pattern="Diagnosis") {
    stopifnot(require(sqldf))
  ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
  pt <- obj$summary$pt
  events <- obj$summary$events
  sqldf("select year, sum(pt) as pt, sum(n) as n, sum(rate*wt) as rate from (select cohort+age as year, age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select cohort, age, sum(pt) as pt from pt group by cohort, age) as t1 natural left outer join (select cohort, age, sum(n) as n from events natural join ev group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year")
}
plotEvents <- function(pattern,ylab="Rate",main=NULL,legend.x="topleft",
                       include.legend=TRUE, legend.y=NULL) {
  with(eventRates(noScreening,pattern),
       plot(year, rate, type="l",ylim=c(0,max(eventRates(goteborg,pattern)$rate)),
            xlab="Age (years)", ylab=ylab, main=main))
  with(eventRates(uptake,pattern), lines(year, rate, col="red"))
  with(eventRates(goteborg,pattern), lines(year, rate, col="green"))
  with(eventRates(riskStrat,pattern), lines(year, rate, col="blue"))
  if (include.legend)
    legend(legend.x,legend.y,
           legend=c("No screening",
             "Opportunistic screening",
             "Göteborg protocol (2+2)",
             "Risk-stratified protocol (4+8)"),
           lty=1,
           col=c("black","red","green","blue"),
           bty="n")
}
prevRatios <- function(obj,predicate) {
  ## ev <- data.frame(event=grep(pattern,levels(obj$summary$prev$event),value=TRUE))
    stopifnot(require(sqldf))
  prev <- obj$summary$prev
  sqldf(sprintf("select year, sum(n) as n, sum(y) as y, sum(p*wt) as prev from (select cohort+age as year, age, t1.n as n, coalesce(t2.y,0.0) as y, 1.0*coalesce(t2.y,0.0)/t1.n*1.0 as p from (select cohort, age, sum(count) as n from prev group by cohort, age) as t1 natural left outer join (select cohort, age, sum(count) as y from prev where %s group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year", predicate))
}
plotPrev <- function(pattern,ylab="Prevalence",main=NULL,legend.x="topleft",
                       include.legend=TRUE, legend.y=NULL) {
  with(prevRatios(noScreening,pattern),
       plot(year, prev, type="l",ylim=c(0,max(prevRatios(goteborg,pattern)$prev)),
            xlab="Age (years)", ylab=ylab, main=main))
  with(prevRatios(uptake,pattern), lines(year, prev, col="red"))
  with(prevRatios(goteborg,pattern), lines(year, prev, col="green"))
  with(prevRatios(riskStrat,pattern), lines(year, prev, col="blue"))
  if (include.legend)
    legend(legend.x,legend.y,
           legend=c("No screening",
             "Opportunistic screening",
             "Göteborg protocol (2+2)",
             "Risk-stratified protocol (4+8)"),
           lty=1,
           col=c("black","red","green","blue"),
           bty="n")
}
table(goteborg$summary$events$event)
table(goteborg$summary$prev$dx)

##path <- function(filename) sprintf("/media/sf_C_DRIVE/usr/tmp/tmp/%s",filename)
##pdf(path("screening_20130425.pdf"),width=7,height=6)
##par(mfrow=c(2,2))
pdf(file="~/work/screening-comparison.pdf",width=6,height=5)
plotEvents("^toScreen$",main="PSA screen",legend.x="topleft")
dev.off()
pdf(file="~/work/biopsy-comparison.pdf",width=6,height=5)
plotEvents("Biopsy",main="Biopsies",legend.x="topleft")
dev.off()
pdf(file="~/work/diagnosis-comparison.pdf",width=6,height=5)
plotEvents("Diagnosis",main="Prostate cancer incidence",legend.x="topleft")
dev.off()
pdf(file="~/work/clinicaldx-comparison.pdf",width=6,height=5)
plotEvents("^toClinicalDiagnosis$",legend.x="bottomleft",
           main="PC incidence (Clinical Dx)")
dev.off()
pdf(file="~/work/prevalence-comparison.pdf",width=6,height=5)
plotPrev("dx!='NotDiagnosed'",main="PC diagnosis",legend.x=2010,legend.y=0.04)
dev.off()
pdf(file="~/work/mortality-comparison.pdf",width=6,height=5)
plotEvents("^toCancerDeath$",legend.x="bottomleft",main="PC mortality")
dev.off()
##dev.off()

## extend the plots to include general conditions
eventRatesCondition <- function(obj,condition,substitute.condition=FALSE) {
  if (substitute.condition)
    condition <- substitute(condition)
  events <- eval(substitute(subset(obj$summary$events, condition), list(condition=condition)))
  pt <- obj$summary$pt
  sqldf("select year, sum(pt) as pt, sum(n) as n, sum(rate*wt) as rate from (select cohort+age as year, age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select cohort, age, sum(pt) as pt from pt group by cohort, age) as t1 natural left outer join (select cohort, age, sum(n) as n from events group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year")
}
plotEventsCondition <- function(condition,ylab="Rate",main=NULL,legend.x="topleft",
                       include.legend=TRUE, legend.y=NULL) {
  condition <- substitute(condition)
  with(eventRatesCondition(noScreening,condition),
       plot(year, rate, type="l",ylim=c(0,max(eventRatesCondition(goteborg,condition)$rate)),
            xlab="Age (years)", ylab=ylab, main=main))
  with(eventRatesCondition(uptake,condition), lines(year, rate, col="red"))
  with(eventRatesCondition(goteborg,condition), lines(year, rate, col="green"))
  with(eventRatesCondition(riskStrat,condition), lines(year, rate, col="blue"))
  if (include.legend)
    legend(legend.x,legend.y,
           legend=c("No screening",
             "Opportunistic screening",
             "Göteborg protocol",
             "Risk-stratified protocol (4+8)"),
           lty=1,
           col=c("black","red","green","blue"),
           bty="n")
}
plotEventsCondition(grepl("Diagnosis",event) & grade %in% c("Gleason_7","Gleason_ge_8"),
                    legend.x="bottomleft",main="PC incidence Gleason 7+")

## How many PSA tests, biopsies etc in the eight years of follow-up?
ratio <- 35000/sum(subset(goteborg$summary$prev,year==2013)$count)
lastYear <- 2014+8
describe <- function(a,b)
  sprintf("pchange=%.1f%%,change=%.1f",100*(1-b/a),(a-b)*ratio)
## PSA tests
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("toScreen",event))$n),
  sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("toScreen",event))$n))
## biopsies
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("toBiopsy",event))$n),
         sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("toBiopsy",event))$n))
## Cancer Gleason 6
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & grade == "Gleason_le_6")$n),
  sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & grade == "Gleason_le_6")$n))
## Cancer Gleason 7+
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & grade %in% c("Gleason_7","Gleason_ge_8"))$n),
  sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & grade %in% c("Gleason_7","Gleason_ge_8"))$n))
## metastatic cancer
describe(sum(subset(riskStrat$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & state=="Metastatic")$n),
  sum(subset(goteborg$summary$events,year>=2013 & year<=lastYear & grepl("Diagnosis",event) & state == "Metastatic")$n))



## In summary, compared with the modified Göteborg protocol over eight years of follow-up, the risk-stratified protocol is expected to have 30% fewer PSA tests (approximately 15,000 fewer), 10% fewer biopsies (~2000 fewer), 3% fewer prostate cancer diagnoses for Gleason 6 cancers (~40 fewer) and 2% fewer cancer diagnoses for Gleason 7+ cancers (~15 fewer cases).
##
## These results are very similar to those obtained by Gulati et al (2013?).
  

  
plotPrev("dx='NotDiagnosed' and state!='Healthy'",main="Latent disease",
         legend.x=2010,legend.y=0.04)
plotPrev("dx='NotDiagnosed' and state!='Healthy' and psa='PSA>=3'",
         main="Latent screen-detectable disease",
         legend.x=2010,legend.y=0.04)

temp0 <- subset(uptake$summary$prev,year==2010 & dx=='NotDiagnosed')
temp <- subset(temp0,state!='Healthy' & psa=='PSA>=3')
temp2 <- sqldf("select pop.age,pop,coalesce(n,0) as n,coalesce(n,0)*1.0/pop as prev  from
(select age, sum(n) as pop from temp0 group by age) as pop natural left join
(select age, sum(n) as n from temp group by age) as cases")
##
temp <- subset(temp0,psa=='PSA>=3')
temp3 <- sqldf("select pop.age,pop,coalesce(n,0) as n,coalesce(n,0)*1.0/pop as prev  from
(select age, sum(n) as pop from temp0 group by age) as pop natural left join
(select age, sum(n) as n from temp group by age) as cases")
with(temp3, plot(age,prev,type="l",ylim=c(0,0.55)))
with(temp2, lines(age,prev,lty=2))

w <- with(subset(pop,age>=50 & age<70),data.frame(age=age,wt=pop/sum(pop)))

sqldf("select sum(wt*prev)/sum(wt) from temp3 natural join w") # prev of PSA 3+ | Not diagnosed
sqldf("select sum(wt*prev)/sum(wt) from temp2 natural join w") # prev of PSA 3+ & cancer | Not diagnosed
5e4*sqldf("select sum(wt*prev)/sum(wt) from temp3 natural join w")
5e4*sqldf("select sum(wt*prev)/sum(wt) from temp2 natural join w")



## Plot of the cohorts over the Lexis diagram
plot(c(1900,2030),c(0,100),type="n",xlab="Calendar period",ylab="Age (years)")
polygon(c(1900,1980,1980+50,1980+50,1980+50-30,1900),
        c(0,0,50,100,100,0))
polygon(c(1990,2030,2030,1990),
        c(50,50,80,80),
        lty=2, border="blue")



### FHCRC model ###
options(width=110)
require(microsimulation)
n <- 1e6
n.cores <- 4
temp <- callFhcrc(n,screen="noScreening",mc.cores=n.cores)
temp4 <- callFhcrc(n,screen="fourYearlyScreen50to70",mc.cores=n.cores)
temp2 <- callFhcrc(n,screen="twoYearlyScreen50to70",mc.cores=n.cores)
temp50 <- callFhcrc(n,screen="screen50",mc.cores=n.cores)
temp60 <- callFhcrc(n,screen="screen60",mc.cores=n.cores)
temp70 <- callFhcrc(n,screen="screen70",mc.cores=n.cores)
uptake <- callFhcrc(n,screen="screenUptake",mc.cores=n.cores)
## "screenUptake", "stockholm3_goteborg", "stockholm3_risk_stratified"

## incremental life-expectancy calculations
LE <- function(obj) sum(obj$summary$pt$pt)/obj$n
IE <- function(obj,objref=temp) LE(obj)-LE(objref)
LE(temp)
IE(temp2)
IE(temp4)
IE(temp50)
IE(temp60)
IE(temp70)

require(data.table)
prev <- data.table(temp2$summary$prev,key="age")
totals <- prev[,sum(count),by="age"]
strat <- prev[,sum(count),by="age,state"]
m <- transform(merge(totals,strat,all=TRUE),prev=V1.y/V1.x)
plot(prev~age+state,m)


require(sqldf)
eventRatesOld <- function(obj,pattern="Diagnosis") {
  ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
  pt <- obj$summary$pt
  events <- obj$summary$events
  sqldf("select year, age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select year, age, sum(pt) as pt from pt group by year, age) as t1 natural left outer join (select year, age, sum(n) as n from events natural join ev group by year, age) as t2")
}
dx <- eventRatesOld(uptake)
require(mgcv)
dx$pred <- 1000*predict(gam(n~s(age,year,k=50),data=dx,offset=log(pt),family=poisson),type="response")
library("RColorBrewer"); library("lattice");
brewer.div <- colorRampPalette(brewer.pal(9, "Spectral"), interpolate = "spline")
pdf("~/work/levelplot-pc-incidence.pdf")
levelplot(pred~age*year,dx,subset=(age>=30), col.regions = brewer.div(100), aspect = "iso",
  xlab="Age (years)", ylab="Calendar period", ylim=c(1980,2050))
dev.off()
with(list(res=600),
       jpeg(file="~/work/levelplot-pc-incidence.jpg",height=5*res,width=5*res,res=res,
            quality=100))
levelplot(pred~age*year,dx,subset=(age>=30), col.regions = brewer.div(100), aspect = "iso",
  xlab="Age (years)", ylab="Calendar period", ylim=c(1980,2050))
dev.off()



## State occupancy: prevalence in different states
prevRatios <- function(obj,predicate) {
  ## ev <- data.frame(event=grep(pattern,levels(obj$summary$prev$event),value=TRUE))
    stopifnot(require(sqldf))
  prev <- obj$summary$prev
  sqldf(sprintf("select year, sum(n) as n, sum(y) as y, sum(p*wt) as prev from (select cohort+age as year, age, t1.n as n, coalesce(t2.y,0.0) as y, 1.0*coalesce(t2.y,0.0)/t1.n*1.0 as p from (select cohort, age, sum(count) as n from prev group by cohort, age) as t1 natural left outer join (select cohort, age, sum(count) as y from prev where %s group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year", predicate))
}
##plotPrev("dx!='NotDiagnosed'",main="PC diagnosis",legend.x=2010,legend.y=0.04)

## rate calculations
## do this using data.table
require(data.table)
eventRates <- function(obj,pattern="Diagnosis") {
  events <- data.table(obj$summary$events,key="event")
  ev <- grep(pattern,levels(events$event),value=TRUE)
  pt <- data.table(obj$summary$pt,key="age")
  with(merge(pt[,sum(pt),by=age],events[J(ev),sum(n),by=age], all=TRUE),
       transform(data.table(age=age,pt=V1.x,n=ifelse(is.na(V1.y),0.0,V1.y))[-1,],
                 rate=n/pt))
}
## the old way: using SQL
require(sqldf)
eventRatesOld <- function(obj,pattern="Diagnosis") {
  ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
  pt <- obj$summary$pt
  events <- obj$summary$events
  sqldf("select age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select age, sum(pt) as pt from pt group by age) as t1 natural left outer join (select age, sum(n) as n from events natural join ev group by age) as t2")
}

## all(abs(eventRates(temp)-data.table(eventRatesOld(temp)))<1e-8)
## system.time(eventRatesOld(temp))
## system.time(eventRates(temp))

png("~/work/screening-comparison-20130222.png",height=2,width=4,res=1200,units="in",pointsize=3)
##x11(width=8,height=5)
##layout(matrix(1:2,nrow=1,byrow=TRUE))
par(mfrow=c(1,2),
  mar      = c(5+2, 4+2, 4+2, 1+2)+0.1,
  ##xaxs     = "i",
  ##yaxs     = "i",
  cex.main = 2,
  cex.axis = 2,
  cex.lab  = 2
)
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Prostate cancer incidence rate",
                            main="None versus two-yearly screening"))
with(eventRates(temp2), lines(age, rate, col="red"))
legend("topleft", legend=c("No screening","Two-yearly\nscreening"), lty=1, col=c("black","red"), bty="n",
       cex=2)
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Prostate cancer incidence rate", main="None versus four-yearly screening"))
with(eventRates(temp4), lines(age, rate, col="blue"))
legend("topleft", legend=c("No screening","Four-yearly\nscreening"), lty=1, col=c("black","blue"), bty="n",
       cex=2)
dev.off()


pdf("~/work/screening-comparison-20130222.pdf")
layout(matrix(1:4,nrow=2,byrow=TRUE))
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus two-yearly screening"))
with(eventRates(temp2), lines(age, rate, col="red"))
legend("topleft", legend=c("No screening","Two-yearly screening"), lty=1, col=c("black","red"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus four-yearly screening"))
with(eventRates(temp4), lines(age, rate, col="blue"))
legend("topleft", legend=c("No screening","Four-yearly screening"), lty=1, col=c("black","blue"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 50"))
with(eventRates(temp50), lines(age, rate, col="green"))
legend("topleft", legend=c("No screening","Screening at age 50"), lty=1, col=c("black","green"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 60"))
with(eventRates(temp60), lines(age, rate, col="orange"))
legend("topleft", legend=c("No screening","Screening at age 60"), lty=1, col=c("black","orange"), bty="n")
dev.off()

pdf("~/work/screening-comparison-20130222.pdf")
layout(matrix(1:4,nrow=2,byrow=TRUE))
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus two-yearly screening"))
with(eventRates(temp2), lines(age, rate, col="red"))
legend("topleft", legend=c("No screening","Two-yearly screening"), lty=1, col=c("black","red"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus four-yearly screening"))
with(eventRates(temp4), lines(age, rate, col="blue"))
legend("topleft", legend=c("No screening","Four-yearly screening"), lty=1, col=c("black","blue"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 50"))
with(eventRates(temp50), lines(age, rate, col="green"))
legend("topleft", legend=c("No screening","Screening at age 50"), lty=1, col=c("black","green"), bty="n")
##
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 60"))
with(eventRates(temp60), lines(age, rate, col="orange"))
legend("topleft", legend=c("No screening","Screening at age 60"), lty=1, col=c("black","orange"), bty="n")
dev.off()

with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None versus screening at age 70"))
with(eventRates(temp70), lines(age, rate, col="orange"))
legend("topleft", legend=c("No screening","Screening at age 70"), lty=1, col=c("black","orange"), bty="n")

personTime <- function(obj) {
  pt <- obj$summary$pt
  n <- obj$n
  sum(pt$pt)/n
}
sapply(list(temp=temp,temp2=temp2,temp4=temp4,temp50=temp50,temp60=temp60,temp70=temp70),personTime) - personTime(temp)
sapply(list(temp=temp,temp2=temp2,temp4=temp4,temp50=temp50,temp60=temp60,temp70=temp70),personTime)




1-exp(-sum(subset(eventRates(temp),age<=75)$rate)) ## 8.4% cumulative risk to age 75 years
1-exp(-sum(subset(eventRates(temp2),age<=75)$rate)) ## 9.4% cumulative risk to age 75 with one random screen between 50 and 70
1-exp(-sum(subset(eventRates(temp3),age<=75)$rate)) ## 9.4% cumulative risk to age 75 with one random screen between 50 and 70

## Stockholm lan, males, 2012

pop <- data.frame(age=0:100,pop=c(12589, 14785, 15373, 14899, 14667,
14437, 14076, 13386, 13425, 12971, 12366, 11659, 11383, 10913, 11059,
11040, 11429, 12303, 13368, 13388, 13670, 13539, 13886, 13913, 14269,
14508, 15073, 15419, 15767, 15721, 16328, 16489, 17126, 16345, 15573,
15702, 16017, 16251, 17069, 16853, 16898, 16506, 15738, 15151, 15224,
15960, 16248, 16272, 16325, 14963, 14091, 13514, 13000, 12758, 12521,
12534, 12333, 11699, 11320, 11167, 11106, 10427, 10889, 10732, 11042,
11367, 11269, 11210, 10982, 10115, 9000, 7652, 6995, 6680, 6144, 5473,
5108, 4721, 4130, 3911, 3756, 3507, 3249, 2803, 2708, 2355, 2188,
2020, 1734, 1558, 1183, 1064, 847, 539, 381, 277, 185, 90, 79, 48,
61))
expected <- function(obj,pattern="Diagnosis",start,end)
sum(transform(subset(merge(obj,eventRates(temp)),age>=50 & age<70),e=rate*pop)$e)



sum(transform(subset(merge(pop,eventRates(temp2)),age>=50 & age<70),e=rate*pop)$e) ## n=748 cases aged 50-69 years (one random screen)
sum(transform(subset(merge(pop,eventRates(temp3)),age>=50 & age<70),e=rate*pop)$e) ## n=748 cases aged 50-69 years (one random screen)
sum(transform(subset(merge(pop,eventRates(temp,"Biopsy")),age>=50 & age<70),e=rate*pop)$e) ## at present, biopsies == cancers
sum(transform(subset(merge(pop,eventRates(temp2,"Biopsy")),age>=50 & age<70),e=rate*pop)$e) ## n=1080 biopsies aged 50-69 years (one random screen)
sum(transform(subset(merge(pop,eventRates(temp3,"Biopsy")),age>=50 & age<70),e=rate*pop)$e) ## n=6920 biopsies aged 50-69 years


plot.EventReport <- function(eventReport, ...) {
  data <- transform(merge(eventReport$pt,eventReport$events), rate=n/pt)
  data <- data[with(data,order(state,event,age)),]
  with(data, plot(age,rate,type="n",ylim=c(0,0.2), ...))
  set <- unique(subset(data,select=c(state,event)))
  invisible(lapply(1:nrow(set),
                   function(i)
                   with(subset(data, state==set$state[i] & event==set$event[i]),
                        lines(age,rate,lty=i))))
} ## Rates only. TODO: add prevalence.
plot.EventReport(temp$summary, xlab="Age (years)")

sapply(temp$parameters,summary)

summary.EventReport <- function(eventReport) {
  data <- transform(merge(eventReport$pt,eventReport$events), rate=n/pt)
  data[with(data,order(state,event,age)),]
}
head(summary.EventReport(temp$summary))


## totals <- data.frame(xtabs(n~age,temp$prev))
## names(totals)[2] <- "total"
## bystate <- data.frame(xtabs(n~age+state,temp$prev))
## names(bystate)[3] <- "n"
## merge(bystate,totals)

require(sqldf)
prev <- temp$prev
temp2 <- sqldf("select *, 1.0*n/total as prev from (select age, sum(n) as total from prev group by age) as t1 left outer join prev on t1.age=prev.age order by state, age")
with(temp2, plot(age,prev,type="n"))
set <- unique(subset(temp2,select=c(state)))
invisible(lapply(1:nrow(set),
                 function(i)
                 with(subset(temp2, state==set$state[i]),
                      lines(age,prev,lty=i))))
legend("bottomleft", legend=set$state, lty=1:7, bty="n")


set.seed(123)
system.time(df <- callSimplePerson(100000))


oldRNGkind <- RNGkind()
if (exists(".Random.seed")) old.Random.seed <- .Random.seed

RNGkind("user")
df <- callSimplePerson()

RNGkind("user")
df <- callPersonSimulation(n=100)

if (exists("old.Random.seed")) .Random.seed <- old.Random.seed
do.call("RNGkind",as.list(oldRNGkind))


require(microsimulation)
Simulation <-
  setRefClass("Simulation",
              contains = "BaseDiscreteEventSimulation",
              fields = list(id = "numeric", state = "character", report = "data.frame"),
              methods= list(initialize = function(id = 0) callSuper(id = id)))
Simulation$methods(init = function() {
  clear()
  id <<- id + 1
  state <<- "Healthy"
  scheduleAt(rweibull(1,8,85), "Death due to other causes")
  scheduleAt(rweibull(1,3,90), "Cancer diagnosis")
})
Simulation$methods(handleMessage = function(event) {
  report <<- rbind(report, data.frame(id = id,
                                      state = state,
                                      begin = previousEventTime,
                                      end = currentTime,
                                      event=event,
                                      stringsAsFactors = FALSE))
  if (event %in% c("Death due to other causes", "Cancer death")) {
    clear()
  }
  else if (event == "Cancer diagnosis") {
    state <<- "Cancer"
    if (runif(1) < 0.5)
      scheduleAt(now() + rweibull(1,2,10), "Cancer death")
  }
})
RNGkind("Mersenne-Twister")
if (exists(".Random.seed")) rm(.Random.seed)
set.seed(123)
sim <- Simulation$new()
system.time(for (i in 1:1000) sim$run())
subset(sim$report,id<=4)

RNGkind("Mersenne-Twister") # cf. "L'Ecuyer-CMRG"!
set.seed(123)
head(.Random.seed)
rng1 <- RNGStream(nextStream = FALSE)
rng2 <- RNGStream()
with(rng1,rexp(1))
with(rng2,rexp(1))
rng1$nextSubStream()
with(rng1,rexp(1))
##
rng1$resetStream()
rng2$resetStream()
with(rng1,rexp(2))
with(rng2,rexp(2))
rng1$nextSubStream()
with(rng1,rexp(2))
rng1$resetRNGkind() # be a good citizen
head(.Random.seed)



temp <- callSimplePerson2(100)
temp2 <- transform(merge(temp$pt,temp$events), rate=n/pt)
temp2 <- temp2[with(temp2,order(state,event,age)),]
with(temp2, plot(age,rate,type="n"))
set <- unique(subset(temp2,select=c(state,event)))
invisible(lapply(1:nrow(set),
                 function(i)
                 with(subset(temp2, state==set$state[i] & event==set$event[i]),
                      lines(age,rate,lty=i))))


muWeibull <- function(a,b) b*gamma(1+1/a)
varWeibull <- function(a,b)  b^2 * (gamma(1 + 2/a) - (gamma(1 + 1/a))^2)
bWeibull <- function(a,mu) mu/gamma(1+1/a)
plotWeibull <- function(a,b,max=60) { x <- seq(0,max,length=301); plot(x,dweibull(x,a,b),type="l") }
##
muWeibull(2,10)
sqrt(varWeibull(2,10))
plotWeibull(2,10)
##
muWeibull(2,3)
sqrt(varWeibull(2,3))
plotWeibull(2,3)
