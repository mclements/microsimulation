## try(detach("package:microsimulation", unload=TRUE))
## require(microsimulation)
## microsimulation:::.testPackage()

## unit tests
refresh
require(microsimulation)
temp=callCalibrationPerson(10)
stopifnot(temp$StateOccupancy[1:2] == c(422,354))
temp2=callFhcrc(1000)
stopifnot(abs(with(temp2,sum(summary$pt$pt)/n)-79.88935)<1e-3)
temp3 <- callIllnessDeath(10)
stopifnot(abs(with(temp3,sum(pt$pt)/10)-64.96217)<1e-3)

temp=callCalibrationPerson(10)
stopifnot(temp$StateOccupancy[1:2] == c(422,354))
temp3 <- callIllnessDeath(10)
stopifnot(abs(with(temp3,sum(pt$pt)/10)-64.96217)<1e-3)

refresh
require(microsimulation)
require(sqldf)
temp2=callFhcrc(1e7,screen="screenUptake",mc.cores=2)

events <- temp2$summary$events
pt <- temp2$summary$pt
m <- dbDriver("SQLite")
connection <- dbConnect(m, dbname = ":memory:")
init_extensions(connection)
dbWriteTable(connection, '`events`', events, row.names = FALSE)
dbWriteTable(connection, '`pt`', pt, row.names = FALSE)

## cancer mortality rates
t1 <- dbGetQuery(connection, "select *, n/pt as rate from (select year, min(85,floor(age/5)*5) as age5, sum(n*1.0) as n from events where event in ('toCancerDeath') and year>=1995 and year<=2010 group by year, age5) t1 natural join (select year, min(85,floor(age/5)*5) as age5, sum(pt) as pt from pt where year>=1995 and year<=2010 group by year, age5) as t2 order by t1.age5, t1.year")

## incidence rates
t1 <- dbGetQuery(connection, "select *, n/pt as rate from (select year, min(85,floor(age/5)*5) as age5, sum(n*1.0) as n from events where event in ('toClinicalDiagnosis','toScreenDiagnosis') and year>=1995 and year<=2010 group by year, age5) t1 natural join (select year, min(85,floor(age/5)*5) as age5, sum(pt) as pt from pt where year>=1995 and year<=2010 group by year, age5) as t2 order by t1.age5, t1.year")
require(lattice)
xyplot(rate ~ year | factor(age5), data=t1, type="l")

dbGetQuery(connection, 'select event,sum(n) from events group by event')

temp=callFhcrc(1e5,nLifeHistories=1e5)
life=temp$lifeHistories
deaths <- sqldf("select t1.id, t1.end as dx, t2.end as death, t1.state, t1.ext_grade from life as t1 inner join life as t2 on t1.id=t2.id where t1.event='toClinicalDiagnosis' and t2.event='toCancerDeath'")

dbDisconnect(connection)

with(subset(deaths,state=="Localised"),plot(density(death-dx,from=0)))

## what proportion of clinical diagnoses die due to cancer?
deaths <- sqldf("select t1.id, t1.end as dx, t2.end as death, t2.event, t1.state, t1.ext_grade from life as t1 inner join life as t2 on t1.id=t2.id where t1.event='toClinicalDiagnosis' and t2.event in ('toOtherDeath','toCancerDeath')")
xtabs(~event+pmin(85,floor(death/5)*5),deaths,subset=(state=="Localised"))

refresh
require(microsimulation)
model <- function(..., mc.cores=3) callFhcrc(..., mc.cores=mc.cores)
model0.0 <- model(1e6,discountRate=0)
model0 <- model(1e6)
model2 <- model(1e6,screen="twoYearlyScreen50to70")
model2.0 <- model(1e6,screen="twoYearlyScreen50to70",discountRate=0)
model4 <- model(1e6,screen="fourYearlyScreen50to70")
ICER(model0.0,model2.0)
ICER(model0,model2)
ICER(model0,model4)
ICER(model2,model4)

## ext_state
refresh
require(microsimulation)
report <- function(obj) {
    rowpct <- function(m) t(apply(m,1,function(row) row/sum(row)))
    print(xtabs(~I(floor(age/5)*5)+ext_state, data=test$diagnoses))
    print(rowpct(xtabs(~I(floor(age/5)*5)+ext_state, data=test$diagnoses)))
    rowpct(xtabs(~I(floor(age/5)*5)+ext_state, data=test$diagnoses, subset=ext_state %in% c('T1_T2','T3plus')))
}
test <- callFhcrc(1e5,includeDiagnoses=TRUE,mc.cores=2)
report(test)
test2 <- callFhcrc(1e5,screen="screenUptake",includeDiagnoses=TRUE,mc.cores=2)
report(test2)
test3 <- callFhcrc(1e5,screen="twoYearlyScreen50to70",includeDiagnoses=TRUE,mc.cores=2)
report(test3)


#### Calibrate for survival
require(microsimulation)
require(dplyr)
## competing risks - using FhcrcParameters$mu0 and fhcrcData$survival_local
CR <- function(agedx,times=c(10,15),HR=1,grade=0) {
    S0 <- data.frame(age=0:105,mu0=FhcrcParameters$mu0) %>%
        filter(age>=agedx) %>%
            mutate(Time=age-agedx,S0=exp(-cumsum(c(0,mu0[-length(mu0)])))) %>%
                select(Time,mu0,S0)
    S1 <- filter(fhcrcData$survival_local,AgeLow==agedx,Grade==grade) %>%
        mutate(S1=Survival^HR,mu1=-HR*log(c(Survival[-1],NA)/Survival)) %>%
            select(Time,mu1,S1)
    inner_join(S0,S1,by="Time") %>%
        mutate(x0=S0*S1*mu0, x1=S0*S1*mu1, CR0=cumsum(x0), CR1=cumsum(x1)) %>%
            filter(Time %in% (times-1)) %>%
                mutate(Time=times) %>%
                    select(Time,CR1,CR0)
}
CRsolve <- function(data) {
    data$grade <- ifelse(data$grade %in% c(0,6,7),0,1)
    optimize(function(hr)
             CR(unique(data$age),HR=hr,grade=data$grade) %>% inner_join(data, by="Time") %>%
             with(sum(c(CR1-cr1)^2)),
             c(0.1,10))
}
## PSA<10, M0/MX
CRsolve(data.frame(age=55,grade=6,Time=c(10,15),cr1=c(0.009,0.029))) # HR=0.137
## CRsolve(data.frame(age=50,grade=7,Time=c(10,15),cr1=c(NA,NA))) # too few at risk
CRsolve(data.frame(age=65,grade=6,Time=c(10,15),cr1=c(0.014,0.047))) # HR=0.181
CRsolve(data.frame(age=65,grade=7,Time=c(10,15),cr1=c(0.087,0.198))) # HR=0.874
CRsolve(data.frame(age=75,grade=6,Time=c(10,15),cr1=c(0.049,0.119))) # HR=0.480
CRsolve(data.frame(age=75,grade=7,Time=c(10,15),cr1=c(0.128,0.217))) # HR=1.020
##
## PSA>=10, M0/MX
CRsolve(data.frame(age=55,grade=6,Time=c(10,15),cr1=c(0.060,0.178))) # HR=0.902
CRsolve(data.frame(age=55,grade=7,Time=c(10,15),cr1=c(0.369,0.474))) # HR=3.631
CRsolve(data.frame(age=65,grade=6,Time=c(10,15),cr1=c(0.088,0.188))) # HR=0.840
CRsolve(data.frame(age=65,grade=7,Time=c(10,15),cr1=c(0.329,0.436))) # HR=2.644
CRsolve(data.frame(age=75,grade=6,Time=c(10,15),cr1=c(0.140,0.239))) # HR=1.134
CRsolve(data.frame(age=75,grade=7,Time=c(10,15),cr1=c(0.265,0.363))) # HR=2.035
CRsolve(data.frame(age=75,grade=8,Time=c(10,15),cr1=c(0.463,0.524))) # HR=0.797
## NB: the SEER survival data were NOT stratified by PSA value.

refresh
require(rstpm2)
require(foreign)
require(lattice)
d <- read.dta("~/Downloads/prostate-20141010.dta")
d2 <- subset(d, age>=60 & age<70 & diayear>=1980)
## debug(pstpm2)
fit <- pstpm2(Surv(time,pcdeath)~1,data=d2,smooth.formula=~s(time,diayear,k=30),sp=0.001)

xtabs(~diayear+addNA(m),data=d)
recent <- subset(d,diayear>=2004 & !is.na(m))
xtabs(~I(floor(age/5)*5)+m,recent)
require(sqldf)
sqldf("select min(90,floor(age/5)*5) as age5, count(*) as n, avg(m) as p_m from recent group by age5")

grid <- expand.grid(diayear=1980:2009,
                    time=seq(0.1,6000,length=50))
grid$haz <- predict(fit,newdata=grid,type="hazard")
grid$surv <- predict(fit,newdata=grid)
xyplot(haz~time | diayear, data=grid, type="l")
xyplot(surv~time | diayear, data=grid, type="l")

xyplot(surv~time | factor(diayear), data=grid,
       panel=function(x,y,subscripts) {
           panel.xyplot(x,y,type="l")
           d3 <- subset(d2,diayear == grid$diayear[subscripts][1])
           sfit <- survfit(Surv(time,pcdeath)~1,data=d3)
           panel.lines(sfit$time,sfit$surv,col=1)
           panel.lines(sfit$time,sfit$lower,col=1,lty=2)
           panel.lines(sfit$time,sfit$upper,col=1,lty=2)
       })

xyplot(pcdeath~time | factor(diayear), data=d,
       subset=age>=50 & age<70,
       panel=function(x,y,subscripts) {
           d3 <- d[subscripts,]
           sfit <- survfit(Surv(time,pcdeath)~1,data=d3)
           panel.lines(sfit$time,sfit$surv,col=1,type="S")
           panel.lines(sfit$time,sfit$lower,col=1,lty=2,type="S")
           panel.lines(sfit$time,sfit$upper,col=1,lty=2,type="S")
       })

xyplot(pcdeath~time | factor(diayear), data=d,
       subset=age>=85 & diayear>=1961,
       xlim=c(0,365*15),
       panel=function(x,y,subscripts) {
           d3 <- d[subscripts,]
           sfit <- survfit(Surv(time,pcdeath)~1,data=d3)
           panel.lines(sfit$time,sfit$surv,col=1,type="S")
           panel.lines(sfit$time,sfit$lower,col=1,lty=2,type="S")
           panel.lines(sfit$time,sfit$upper,col=1,lty=2,type="S")
       })




## all(c(callFhcrc(5)$lifeHistories == callFhcrc(5)$lifeHistories,
##       callFhcrc(5)$parameters == callFhcrc(5)$parameters,
##       callFhcrc(5)$summary$pt == callFhcrc(5)$summary$pt,
##       callFhcrc(5)$summary$prev == callFhcrc(5)$summary$prev))
## all(callFhcrc(10)$lifeHistories == callFhcrc(10)$lifeHistories) # fails
## all(callFhcrc(1e4)$parameters == callFhcrc(1e4)$parameters) # okay

## extract PSA values and calculate STHLM3 PSA "pseudo-thresholds" for the biomarker panel
refresh
require(microsimulation)
pos <- function(x) ifelse(x>0,x,0)
set.seed(12345)
temp <- callFhcrc(1e6,screen="stockholm3_risk_stratified",includePSArecords=TRUE,mc.cores=2)$psarecord
temp2 <- subset(temp,organised & age>=50 & age<70 & !dx) 
## first organised screen
i <- tapply(1:nrow(temp2),temp2$id,min)
temp3 <- temp2[i,]
cat("No cancer:\n")
with(subset(temp3,state==0 & psa>3), mean(psa<4.4)) # about 42% (from STHLM3)
cat("Loco-regional Cancer:\n")
with(subset(temp3,state>0 & ext_grade==0 & psa>3), mean(psa<3.6)) # about 17% (from STHLM3)
with(subset(temp3,state>0 & ext_grade==0 & psa>3),
     cumsum(table(cut(delta,c(0,5,10,15,Inf)))/length(delta)))
with(subset(temp3,state>0 & ext_grade==0 & psa>3.6), 
     plot(density(delta,from=0)))
with(subset(temp3,state>0 & ext_grade==0 & psa>3), 
    lines(density(delta,from=0),lty=2))

with(subset(temp3,state>0 & ext_grade==0 & psa>3), mean(delta))
with(subset(temp3,state>0 & ext_grade==0 & psa>3 & psa<3.6), mean(delta))
temp3 <- transform(temp3, delta=age-(t0+35))


temp2 <- transform(temp2,
                   advanced=(state>0 & ext_grade==2),
                   cancer=(state>0),
                   logpsa=log(psa),
                   logZ=log(Z),
                   logZstar=beta0+beta1*(age-35)+pos(beta2*(age-35-t0)),
                   grade=ifelse(ext_grade %in% 0:1,1,ext_grade))
temp2 <- within(temp2, {
    ext_grade <- ifelse(cancer, ext_grade, NA)
})
## ## rnormPos is now in the package...
## rnormPos <- function(n,mean=0,sd=1,lbound=0) {
##     if (length(mean)<n) mean <- rep(mean,length=n)
##     if (length(sd)<n) sd <- rep(sd,length=n)
##     x <- rnorm(n,mean,sd)
##     while(any(i <- which(x<lbound)))
##         x[i] <- rnorm(length(i),mean[i],sd[i])
##     x
## }
## onset ho(t) = g0 * t
p <- list(mubeta0=-1.609,
          sebeta0=0.2384,
          mubeta1=0.04463,
          sebeta1=0.0430,
          mubeta2=c(0.0397,0.1678),
          sebeta2=c(0.0913,0.3968),
          tau2=0.0829,
          g0=0.0005)
n <- nrow(temp2)
correlatedValue <- function(y,mu,se=NULL,rho) {
    residual <- y - mu 
    if (is.null(se)) se <- sqrt(var(residual))
    u <- rnorm(length(y),0,se) # marginal error
    mu + rho*residual + sqrt(1-rho^2)*u # new value
}

rtpf <- function(marker1, threshold1, marker2, threshold2, disease) {
    n1 <- sum(marker1[disease] > threshold1)
    n2 <- sum(marker2[disease] > threshold2)
    list(n1=n1, n2=n2, rtpf=n1/n2)
}
rfpf <- function(marker1, threshold1, marker2, threshold2, disease) {
    n1 <- sum(marker1[!disease] > threshold1)
    n2 <- sum(marker2[!disease] > threshold2)
    list(n1=n1, n2=n2, rfpf=n1/n2)
}
variance <- tau2 <- 0.0829
optim1 <- optim(c(log(4.5),log(0.01)),
                function(par) {
                    set.seed(12345)
                    threshold <- par[1]
                    scale <- exp(par[2])
                    temp2$bbp <<- with(temp2, logpsa + scale*pos(age-t0) + rnorm(nrow(temp2), 0, sqrt(variance)))
                    with(temp2,
                         (rtpf(logpsa, log(3.0), bbp, threshold, advanced)$rtpf-1)^2 +
                         (rfpf(logpsa, log(3.0), bbp, threshold, advanced)$rfpf-1.25)^2)
                })
set.seed(12345)
threshold <- optim1$par[1]
scale <- exp(optim1$par[2])
temp2$bbp <- with(temp2, logpsa + scale*pos(age-t0) + rnorm(nrow(temp2), 0, sqrt(variance)))
with(temp2, rtpf(logpsa, log(3.0), bbp, threshold, advanced))
with(temp2, rfpf(logpsa, log(3.0), bbp, threshold, advanced))

with(temp2, rtpf(logpsa, log(3.0), bbp, log(10), advanced))
with(temp2, rfpf(logpsa, log(3.0), bbp, log(10), advanced))
root1 <- uniroot(function(threshold) with(temp2, rtpf(logpsa, log(3.0), bbp, threshold, advanced))$rtpf-1, c(log(2), log(100)))
with(temp2, rtpf(logpsa, log(3.0), bbp, root1$root, advanced))
with(temp2, rfpf(logpsa, log(3.0), bbp, root1$root, advanced))


## correlated biomarkers based on the mean
set.seed(12345+5)
biomarker2 <- exp(correlatedValue(log(temp2$psa),
                                  p$mubeta0+p$mubeta1*(temp2$age-35)+p$mubeta2[temp2$grade]*pos(temp2$age-35-temp2$t0),
                                  rho=0.25))
biomarker2 <- exp(log(temp2$psa) + 0.1*pos(temp2$age-35-temp2$t0)+rnorm(n,0,sqrt(p$tau2)))
if (FALSE) {
    plot(temp2$psa,biomarker2,log="xy")
    sqrt(var(log(temp2$psa) - p$mubeta0+p$mubeta1*(temp2$age-35)+p$mubeta2*pos(temp2$age-35-temp2$t0)))
    cor(log(biomarker2),log(temp2$psa))
    plot(density(log(temp2$psa)))
    lines(density(log(biomarker2)),lty=2)
    var(log(biomarker2))
    var(log(temp2$psa))
}
temp3 <- transform(temp2, BBP=biomarker2)
## STHLM3 simulation report
with(list(threshold=5),with(transform(temp3,BBPpos=(psa>=1 & BBP>=threshold),PSApos=(psa>=3)),
     cat(sprintf("
PSA+ & advanced:\t\t%i
PSA+ & cancer:\t\t\t%i
BBP+ & adv:\t\t\t%i
BBP+ & can:\t\t\t%i
BBP+ & PSA+:\t\t\t%i
BBP+ & PSA-:\t\t\t%i
BBP- & PSA+:\t\t\t%i
PSA+ | BBP+:\t\t\t%i
Prop reduction in biospies:\t%5.3f\n",
                 sum(PSApos & advanced),
                 sum(PSApos & cancer),
                 sum(BBPpos & advanced),
                 sum(BBPpos & cancer),
                 sum(BBPpos & PSApos),
                 sum(BBPpos & !PSApos),
                 sum(!BBPpos & PSApos),
                 sum(BBPpos | PSApos),
                 (sum(!BBPpos & PSApos) - sum(BBPpos & !PSApos))/
                       sum(PSApos)
                 ))))

report <- function(psa, BBP, advanced, threshold=3) {
    BBPpos <- (psa>=1 & BBP>=threshold)
    PSApos <- (psa>=3)
    c(PSAposAdvanced=sum(PSApos & advanced),
      BBPposAdvanced=sum(BBPpos & advanced),
      pBiopsy=(sum(!BBPpos & PSApos) - sum(BBPpos & !PSApos))/
      sum(PSApos))
}
report(temp2$psa,biomarker2,temp2$advanced,threshold=1.8)
reports <- sapply(1:200,function(i) {
    biomarker2 <- exp(correlatedValue(log(temp2$psa),
                                      p$mubeta0+p$mubeta1*(temp2$age-35)+p$mubeta2[temp2$grade]*pos(temp2$age-35-temp2$t0),
                                      rho=0.25))
    report(temp2$psa,biomarker2,temp2$advanced, threshold=1.8)
})
reports <- as.data.frame(t(reports))
with(reports, mean( PSAposAdvanced <= BBPposAdvanced))
with(reports, plot(table( PSAposAdvanced - BBPposAdvanced)))
with(reports, plot(density(pBiopsy)))
with(reports, plot(PSAposAdvanced - BBPposAdvanced,pBiopsy))

## Old natural history - ignoring diagnosis
p <- list(mubeta0=-1.609,
          sebeta0=0.2384,
          mubeta1=0.04463,
          sebeta1=0.0430,
          mubeta2=c(0.0397,0.1678),
          sebeta2=c(0.0913,0.3968),
          tau2=0.0829,
          g0=0.0005)
## onset ho(t) = g0 * t, Ho(t) = g0/2*t*t = -logU => t=sqrt(-2*log(U)/g0)
set.seed(12345)
n <- 1e6
age_o <- 35+sqrt(-2*log(runif(n))/p$g0)
## grade <- rep(1:2,c(0.9*n,0.1*n)) # this should depend on age of onset
grade <- ifelse(runif(n) < 0.006*(age_o-35), 1, 0)
beta0 <- with(p, rnorm(n,mubeta0,sebeta0))
beta1 <- with(p, rnormPos(n,mubeta1,sebeta1))
beta2 <- with(p, rnormPos(n,mubeta2[grade],sebeta2[grade]))
eps <- with(p, rnorm(n,0,tau2))
lpsa <- pmin(log(20),beta0+beta1*(50-35)+beta2*pmax(0,50-age_o)+eps)
psacut <- function(x) cut(x,c(0,1,3,10,Inf), right=FALSE)
table(psacut(exp(lpsa)))/length(lpsa)
plot(density(exp(lpsa)),xlim=c(0,20)) # density of PSA at age 50 years
i <- 1
for (age in seq(45,85,by=10)) {
    cat(age,"\n")
    lpsa <- pmin(log(20),beta0+beta1*(age-35)+beta2*pmax(0,age-age_o)+eps)
    lines(density(exp(lpsa)),col=i)
    print(table(cut(exp(lpsa),psa_cuts))/length(lpsa))
    i <- i+1
}
tab <- sapply(ages <- seq(55,80,by=5), function(age) {
    lpsa <- pmin(log(20),beta0+beta1*(age-35)+beta2*pmax(0,age-age_o)+eps)
    lpsa <- pmin(log(20),beta0+beta1*(age-35)+beta2*pmax(0,age-age_o))
    tab <- table(psacut(exp(lpsa)))
    tab/sum(tab)
})
colnames(tab) <- ages
tab



## revised natural history?

p <- list(mubeta0=-1.609,
          sebeta0=0.2384,
          mubeta1=0.04463,
          sebeta1=0.0430,
          mubeta2=c(0.0397,0.1678),
          sebeta2=c(0.0913,0.3968),
          tau2=0.0829,
          g0=0.0005)
psaSim <- function(n,age=50,grade=NULL) {
    age_o <- 35+sqrt(-2*log(runif(n))/p$g0)
    if (is.null(grade)) grade <- rep(1:2,c(0.9*n,0.1*n))
    grade <- rep(grade,length=n)
    beta0 <- with(p, rnorm(n,mubeta0,sebeta0))
    beta1 <- with(p, rnormPos(n,mubeta1,sebeta1))
    beta2 <- with(p, rnormPos(n,mubeta2[grade],sebeta2[grade]))
    eps <- with(p, rnorm(n,0,tau2))
    lpsa <- pmin(log(20),beta0+beta1*(age-35)+beta2*pmax(0,age-age_o)+eps)
    psa <- exp(lpsa)
    as.data.frame(as.list(environment()))
}
set.seed(12345)
psaSim(11)

refresh
require(microsimulation)
y <- t(replicate(1000,.Call("rbinormPos_test",package="microsimulation")))
cor(y)
plot(y)


refresh
require(microsimulation)
require(sqldf)
require(dplyr)
load("~/Downloads/IHEdata.RData")
makeModel <- function(discountRate=0.03,
                      formal_compliance=1,
                      formal_costs=1,
                      panel=TRUE) {
    function(screen, ..., parms=NULL, n=1e6, mc.cores=3, pop=1960) {
        newparms <- list(formal_compliance=formal_compliance,
                         formal_costs=formal_costs)
        if (!is.null(parms))
            for (name in names(parms))
                newparms[[name]] <- parms[[name]]
        callFhcrc(n, screen=screen, mc.cores=mc.cores, pop=pop, discountRate=discountRate, parms=newparms, panel=panel, ...)
    }
}
modelSet <- function(model) {
    model0 <- model("noScreening")
    model2 <- model("twoYearlyScreen50to70")
    model4 <- model("fourYearlyScreen50to70")
    model50 <- model("screen50")
    model60 <- model("screen60")
    model70 <- model("screen70")
    modelUptake1930 <- model("screenUptake",pop=1930)
    modelUptake1960 <- model("screenUptake")
    modelGoteborg <- model("goteborg")
    modelRiskStratified <- model("risk_stratified")
    modelMixedScreening <- model("mixed_screening")
    cat("NOTE: Processing completed.\n")
    models <- list(model0,model2,model4,model50,model60,model70,
                   modelUptake1930,modelUptake1960,modelGoteborg,
                   modelRiskStratified,modelMixedScreening)
    names(models) <- c("No screening","2-yearly","4-yearly",
                       "50 only","60 only","70 only","Opportunistic 1930",
                       "Opportunistic 1960+",
                       "Göteborg","Risk stratified",
                         "Mixed screening")
    models
}
predict.fhcrc <-
function (obj, type = c("incidence", "cancerdeath"))
{
    type <- match.arg(type)
    event_types <- switch(type, incidence = c("toClinicalDiagnosis", 
        "toScreenDiagnosis"), cancerdeath = "toCancerDeath")
    if (require(dplyr)) {
        pt <- obj$summary$pt %>% group_by(age) %>% summarise(pt = sum(pt))
        events <- obj$summary$events %>% filter(event %in% event_types) %>% 
            group_by(age) %>% summarise(n = sum(n))
        out <- left_join(pt, events, by = "age") %>% mutate(rate = ifelse(is.na(n), 
            0, n/pt))
        return(with(out,data.frame(age=age,rate=rate,pt=pt,n=ifelse(is.na(n),0,n))))
    }
    else error("dplyr is not available for plotting")
}
plot.scenarios <- function(models,
                           costs="delta.costs",
                           effects="delta.QALE",
                           xlim=NULL, ylim=NULL,
                           ylab="Effectiveness (QALY)",
                           suffix="", prefix="",
                           textp=TRUE,
                           pos=rep(4,length(models)),
                           ...) {
    s <- data.frame(t(sapply(models,
                             function(obj) unlist(ICER(obj,models[[1]])))),
                    model=sprintf("%s%s%s",prefix,names(models),suffix),
                    pos=pos)
    costs <- s[[costs]]
    effects <- s[[effects]]
    plot(costs,
         effects,
         xlim=if (is.null(xlim)) c(0,max(costs)*1.3) else xlim,
         ylim=if (is.null(ylim)) c(0,max(effects)*1.1) else ylim,
         xlab="Costs (SEK)",
         ylab=ylab,
         pch=19, cex=1.5,
         ...)
    if (textp) text(costs,effects, labels=s$model, pos=pos)
    lines.frontier(costs,effects,type="c",lwd=2,col="grey")
}
points.scenarios <- function(models,
                             costs="delta.costs",
                             effects="delta.QALE",
                             suffix="",
                             prefix="",
                             textp = TRUE,
                             pos=rep(4,length(models)), ...) {
    s <- data.frame(t(sapply(models,
                             function(obj) unlist(ICER(obj,models[[1]])))),
                    model=sprintf("%s%s%s",prefix,names(models),suffix),
                    pos=pos)
    costs <- s[[costs]]
    effects <- s[[effects]]
    points(costs,
           effects,
           pch=19,cex=1.5,
           ...)
    if (textp) text(costs,effects, labels=s$model, pos=pos)
}
segments.scenarios <- function(modelsA,
                               modelsB,
                               costs="delta.costs",
                               textp=FALSE,
                               pos=rep(4,length(modelsA)),
                               effects="delta.QALE",
                               ...) {
    sA <- data.frame(t(sapply(modelsA,
                              function(obj) unlist(ICER(obj,modelsA[[1]])))))
    sB <- data.frame(t(sapply(modelsB,
                              function(obj) unlist(ICER(obj,modelsB[[1]])))))
    costsA <- sA[[costs]]
    effectsA <- sA[[effects]]
    costsB <- sB[[costs]]
    effectsB <- sB[[effects]]
    segments(costsA,effectsA,
             costsB,effectsB,
             lwd=2,
             ...)
    if (textp)
        text((costsA+costsB)/2,
             (effectsA+effectsB)/2,
             labels=names(modelsA),
             pos=pos)
}
summary.scenarios <- function(models) {
    data.frame(t(sapply(models,
                             function(obj) unlist(ICER(obj,models[[1]])))),
                    model=names(models))
}
post <- function(modelSet) {
    i <- c(1,4,5,6,8:11)
    names(modelSet) <- c("No screening","Göteborg","4-yearly",
                         "50 only","60 only","70 only","Opportunistic 1930",
                         "Opportunistic",
                         "Risk stratified (2+4)","Risk stratified (4+8)",
                         "Mixed screening")
    modelSet[i]
}
if (FALSE) {
    modelSetA <- modelSet(makeModel(discount=0,formal_compliance=0,formal_costs=0,panel=FALSE))
    modelSetB <- modelSet(makeModel(discount=0.03,formal_compliance=0,formal_costs=0,panel=FALSE))
    modelSetC <- modelSet(makeModel(discount=0.03,formal_compliance=1,formal_costs=1,panel=FALSE))
    modelSetD <- modelSet(makeModel(discount=0.03,formal_compliance=1,formal_costs=1,panel=TRUE))
    modelSetBD <- modelSet(makeModel(discount=0.03,formal_compliance=0,formal_costs=0,panel=TRUE))
    save(modelSetA,file="~/work/modelSetA-20150201.RData")
    save(modelSetB,file="~/work/modelSetB-20150201.RData")
    save(modelSetC,file="~/work/modelSetC-20150201.RData")
    save(modelSetD,file="~/work/modelSetD-20150201.RData")
    save(modelSetBD,file="~/work/modelSetBD-20150201.RData")
}
doOnce <- TRUE
if (doOnce) {
    load("~/work/modelSetA-20150201.RData")
    load("~/work/modelSetB-20150201.RData")
    load("~/work/modelSetC-20150201.RData")
    load("~/work/modelSetD-20150201.RData")
    load("~/work/modelSetBD-20150201.RData")
    modelSetA <- post(modelSetA)
    modelSetB <- post(modelSetB)
    modelSetC <- post(modelSetC)
    modelSetD <- post(modelSetD)
    modelSetBD <- post(modelSetBD)
    doOnce <- FALSE
}
## ## labels
##     c("1"="No screening",
##       "2"="50 only",
##       "3"="60 only",
##       "4"="70 only",
##       "5"="Opportunistic",
##       "6"="Risk stratified (2+4)",
##       "7"="Risk stratified (4+8)",
##       "8"="Mixed screening")

## modelSetBD0 <- modelSet(makeModel(discount=0,formal_compliance=0,formal_costs=0,panel=TRUE))
## modelSetBD0 <- post(modelSetBD0)
plot.scenarios(c(modelSetA,modelSetBD0),
               type="n",textp=FALSE)
points.scenarios(modelSetBD0,col="violet",textp=FALSE)
points.scenarios(modelSetA,col="red",textp=FALSE)
segments.scenarios(modelSetBD0, modelSetA,textp=TRUE,pos=c(4,1,4,4,1,3,2,1))
legend("bottomright",legend=c("Panel + informal","PSA + informal"),col=c("violet","red"),bty="n",pch=19,pt.cex=1.5)

## Tables of costs and effectiveness
rbind(transform(summary.scenarios(modelSetB),set="B"),
      transform(summary.scenarios(modelSetC),set="C"),
      transform(summary.scenarios(modelSetD),set="D"))

## individual plots
pdf("~/Downloads/cea-A-LY.pdf")
plot.scenarios(modelSetA,effects="delta.LE",ylab="Effectiveness (LY)")
dev.off()
pdf("~/Downloads/cea-A-QALY.pdf")
plot.scenarios(modelSetA)
dev.off()
pdf("~/Downloads/cea-B.pdf")
plot.scenarios(modelSetB,col="red")
dev.off()
pdf("~/Downloads/cea-C.pdf")
plot.scenarios(modelSetC,col="orange")
dev.off()
pdf("~/Downloads/cea-D.pdf")
plot.scenarios(modelSetD,col="green")
dev.off()

## sets B, BD, C and D together
plot.scenarios(c(modelSetB,modelSetC,modelSetD,modelSetBD),type="n",textp=FALSE)
points.scenarios(modelSetC,col="orange")
points.scenarios(modelSetB,col="red")
points.scenarios(modelSetD,col="green",textp=FALSE)
points.scenarios(modelSetBD,col="violet",textp=FALSE)
legend("bottomright",legend=c("Panel + formal","PSA + formal","Panel + informal","PSA + informal"),col=c("green","orange","violet","red"),bty="n",pch=19,pt.cex=1.5)

pdf("~/Downloads/cea-BC.pdf")
plot.scenarios(c(modelSetC,modelSetB),xlim=c(0,3000),
               type="n",textp=FALSE)
points.scenarios(modelSetC,col="orange",textp=FALSE)
points.scenarios(modelSetB,col="red",textp=FALSE)
segments.scenarios(modelSetB, modelSetC,textp=TRUE,pos=c(4,1,4,4,1,3,2,1))
legend("bottomright",legend=c("PSA + formal","PSA + informal"),col=c("orange","red"),bty="n",pch=19,pt.cex=1.5)
dev.off()

pdf("~/Downloads/cea-BvCD.pdf")
plot.scenarios(c(modelSetBD,modelSetB),xlim=c(0,3000),
               type="n",textp=FALSE)
points.scenarios(modelSetBD,col="violet",textp=FALSE)
points.scenarios(modelSetB,col="red",textp=FALSE)
segments.scenarios(modelSetB, modelSetBD,textp=TRUE,pos=c(4,1,4,4,1,3,2,1))
legend("bottomright",legend=c("Panel + informal","PSA + informal"),col=c("violet","red"),bty="n",pch=19,pt.cex=1.5)
dev.off()

pdf("~/Downloads/cea-CD.pdf")
plot.scenarios(c(modelSetC,modelSetD),xlim=c(0,3000),col="orange",textp=FALSE,type="n")
points.scenarios(modelSetC,col="orange",textp=FALSE)
points.scenarios(modelSetD,col="green",textp=FALSE)
segments.scenarios(modelSetC, modelSetD,textp=TRUE)
legend("bottomright",legend=c("PSA + formal","Panel + formal"),col=c("orange","green"),bty="n",pch=19,pt.cex=1.5)
dev.off()

pdf("~/Downloads/cea-BDvD.pdf")
plot.scenarios(c(modelSetD,modelSetBD),xlim=c(0,3000),textp=FALSE,type="n")
points.scenarios(modelSetD,col="green",textp=FALSE)
points.scenarios(modelSetBD,col="violet",textp=FALSE)
segments.scenarios(modelSetBD, modelSetD,textp=TRUE)
legend("bottomright",legend=c("Panel + informal","Panel + formal"),col=c("violet","green"),bty="n",pch=19,pt.cex=1.5)
dev.off()

pdf("~/Downloads/cea-incidence.pdf")
plot(modelSetA[[1]],type="incidence",xlab="Age (years)", ylab="Incidence rate", lwd=2)
for (i in 2:length(modelSetA))
    lines(modelSetA[[i]],col=i,type="incidence", lwd=2)
legend("bottomright",legend=names(modelSetA),lty=1,col=1:length(modelSetA),bty="n",lwd=2)
dev.off()


## Just compliance
model <- makeModel(discount=0,formal_compliance=1,formal_costs=0,panel=FALSE)
model0 <- model("noScreening")
model2 <- model("twoYearlyScreen50to70")
plot(model0,type="cancerdeath")
comparison1 <- rbind(predict(model2,type="cancerdeath") %>% mutate(screen=1),
                    predict(model0,type="cancerdeath") %>% mutate(screen=0)) %>%
    filter(age >= 50 & age<70)
require(splines)
exp(coef(glm(n ~ offset(log(pt)) + ns(age) + screen, data=comparison1, family=poisson)))
##
comparison <-
    inner_join(predict(model2,type="cancerdeath") %>% rename(rate2=rate) %>% select(age,rate2),
               predict(model0,type="cancerdeath") %>% rename(rate0=rate) %>% select(age,rate0)) %>%
    mutate(RR=rate2/rate0) %>% filter(age>=50 & age<90)
plot(RR ~ age, data=comparison, type="l")

## Treatment patterns
treat <- with(stockholmTreatment,
              data.frame(data.frame(year=DxY,Age=Age,G=factor(G+1,labels=c("Gleason 6","Gleason 7","Gleason 8+"))),
                         Treatment=factor(rep(1:3,each=nrow(stockholmTreatment)),labels=c("CM","RP","RT")),
                         Proportion=as.vector(cbind(CM,RP,RT))))
require(ggplot2)
pdf("~/work/treatment_patterns.pdf")
ggplot(treat, aes(x=Age,y=Proportion,group=Treatment,fill=Treatment)) + facet_wrap(~G) + geom_area(position="fill") + xlab("Age (years)")
dev.off()

if (FALSE) {
    plot(modelSetD[["No screening"]],type="incidence")
    lines(modelSetD[[""]],type="incidence",col="blue")
    lines(modelMixedScreening,type="incidence",col="red")
    lines(modelGoteborg,type="incidence",col="green")
    lines(modelRiskStratified,type="incidence",col="lightblue")
    lines(model1,type="incidence",col="orange")
    lines(modelUptake1960,type="incidence",col="pink")
}

## List of homogeneous elements
List <- function(...) {
    .Data <- list(...)
    class.element <- class(.Data[[1]])
    stopifnot(all(sapply(.Data, function(element) class(element)==class.element)))
    structure(.Data=.Data,
              element.class=class.element, # new attribute
              names=names(.Data),
              class = c("List","list"))
}
print.List <- function(obj,...) {
    i <- 1
    namess <- names(obj)
    for (obji in obj) {
        name <- if (is.null(namess)) sprintf("[[%i]]",i) else namess[i]
        cat(name,"\n")
        print(obji,...)
        i <- i+1
    }
    invisible(obj)
}
getListFunction <- function(fun,obj,...) {
    stopifnot(inherits(obj,"List"))
    VALUE <- sprintf("%s.List.%s",
                     deparse(substitute(fun)),
                     attr(obj,"element.class"))
    FUN <- tryCatch(get(VALUE,...))
    if (inherits(FUN,"try-error")) stop(sprintf("%s is not defined.\n",VALUE))
    FUN
}
plot.List <- function(obj,...) {
    getListFunction(plot,obj)(obj,...)
}
plot.List.fhcrc <- function(obj,...) {
    temp <- data.frame(t(sapply(obj,
                                function(obji) unlist(ICER(obji,obj[[1]])))))
    with(temp,
         plot(delta.costs,
              delta.QALE,
              xlim=c(0,max(delta.costs)*1.3),
              ylim=c(0,max(delta.QALE)*1.1),
              xlab="Costs (SEK)",
              ylab="Effectiveness (QALY)"))
    lines.frontier(temp$delta.costs,temp$delta.QALE)
    invisible(obj)
}
plot(List(model0,model1,model2,model50,model60,model70,
          modelUptake1930,modelUptake1960,modelGoteborg,
          modelRiskStratified,modelMixedScreening,modelFormalTestManagement))
List(model0,model1,model2,model50,model60,model70,
     modelUptake1930,modelUptake1960,modelGoteborg,
     modelRiskStratified,modelMixedScreening,modelFormalTestManagement)



## save(model0,model1,model1p,model2,model2p,model50,model60,model70,
##      modelUptake1930,modelUptake1960,modelGoteborg,
##      modelRiskStratified,modelMixedScreening,modelFormalTestManagement,
##      file="~/work/icer_20150111.RData")


mubeta2 <- c(0.0397,0.1678)
p <- 0.9 
fun <- function(par,a,b) {
    alpha <- par[1]
    beta <- par[2]
    (exp(alpha+beta*b)-exp(alpha+beta*a))/(b-a)/beta
}
objective <- function(par) {
    alpha <- par[1]
    beta <- par[2]
    (mubeta2[1]-fun(par,0,p))^2+
    (mubeta2[2]-fun(par,p,1))^2
}
fit <- optim(c(1,1),objective,control=list(abstol=1e-16,reltol=1e-16))
fun(fit$par, p, 1)
fun(fit$par, 0, p)
p1 <- 0.6
fun(fit$par, 0, p1)
fun(fit$par, p1, p)
with(list(x=seq(0,1,length=301)),
     plot(x,sapply(x,function(xi) exp(fit$par[1]+fit$par[2]*xi)), type="l"))
abline(v=p,lty=2)
abline(v=p1, lty=2)
segments(0,mubeta2[1],p,mubeta2[1],lty=3)
segments(p,mubeta2[2],1,mubeta2[2],lty=3)

## Even width
mubeta2 <- c(0.0397,0.1678)
fun <- function(par,a,b) {
    alpha <- par[1]
    beta <- par[2]
    (exp(alpha+beta*b)-exp(alpha+beta*a))/(b-a)/beta
}
objective <- function(par) {
    alpha <- par[1]
    beta <- par[2]
    (mubeta2[1]-fun(par,6,8))^2+
    (mubeta2[2]-fun(par,8,9))^2
}
fit <- optim(c(1,1),objective,control=list(abstol=1e-16,reltol=1e-16))
fun(fit$par, 6, 8)
fun(fit$par, 8, 9)
fun(fit$par, 6, 7)
fun(fit$par, 7, 8)
with(list(x=seq(6,9,length=301)),
     plot(x,sapply(x,function(xi) exp(fit$par[1]+fit$par[2]*xi)), type="l"))


## check cost calculations
model0 <- callFhcrc(1e5,screen="twoYearlyScreen50to70",mc.cores=3,pop=1995-50,discountRate=0)
model1 <- callFhcrc(1e5,screen="fourYearlyScreen50to70",mc.cores=3,pop=1995-50,discountRate=0)
costs <- model1$costs
pt <- model1$summary$pt
pop1 <- sqldf("select age, sum(pt) as pop from pt group by age")
costs1 <- sqldf("select age, item, sum(costs) as costs from costs group by age, item")
sqldf("select item, sum(costs/pop*IHE/1e6) as adj from pop1 natural join costs1 natural join IHEpop group by item")


## Correlated PSA values
refresh
require(microsimulation)
require(mvtnorm)
require(dplyr)
p <- list(mubeta0=-1.609,
          sebeta0=0.2384,
          mubeta1=0.04463,
          sebeta1=0.0430,
          mubeta2=c(0.0397,0.1678),
          sebeta2=c(0.0913,0.3968),
          tau2=0.0829,
          g0=0.0005)
rmvnormPos <- function(n,mean=0,sigma=matrix(1),lbound=0) {
    x <- rmvnorm(n,mean,sigma)
    while(any(i <- which(apply(x,1,min) < lbound)))
        x[i,] <- rmvnorm(length(i),mean,sigma)
    x
}
## rmvnormPos(10,c(0,0),matrix(c(1,0,0,1),2))
prob_grade7 <- fhcrcData$prob_grade7 %>% "names<-"(c("x","y")) %>% approxfun()
psaSimCor <- function(n,age=50,rho=0.62,max.psa=50,mubeta2.scale) {
    age_o <- 35+sqrt(-2*log(runif(n))/p$g0)
    grade <- ifelse(runif(n)>=1+FhcrcParameters$c_low_grade_slope*(age_o-35),8,7)
    cor0 <- cor1 <- cor2 <- matrix(c(1,rho,rho,1),2)
    beta0 <- with(p, rmvnorm(n,c(mubeta0,mubeta0),sebeta0^2*cor0))
    beta1 <- with(p, rmvnorm(n,c(mubeta1,mubeta1),sebeta1^2*cor1))
    beta2 <- matrix(NA,n,2)
    for (gradei in 7:8) {
        i <- which(grade == gradei)
        index <- if(gradei==7) 1 else 2
        if (any(i)) {
            x <- with(p, rmvnorm(length(i),c(mubeta2[index],mubeta2.scale*mubeta2[index]),sebeta2[index]^2*cor2))
            beta2[i,1] <- x[,1]
            beta2[i,2] <- x[,2]
        }
    }
    ext_grade <- ifelse(grade==7,
                        ifelse(runif(n)<prob_grade7(beta2),7,6),
                        8)
    eps <- with(p, cbind(rnorm(n,0,tau2),rnorm(n,0,tau2)))
    lpsa <- t(apply(beta0+beta1*(age-35)+beta2*pmax(0,age-age_o)+eps, 1, pmin, log(max.psa)))
    ## psa <- exp(lpsa)
    data.frame(age=age,
               cancer=age_o<age,
               advCancer=age_o<age, ### !!!!!!!
               lpsa=lpsa[,1],
               lbp=lpsa[,2],
               psa=exp(lpsa[,1]),
               bp=exp(lpsa[,2]),
               age_o=age_o,
               grade=grade,
               ext_grade=ext_grade)
}
dAgg <- function(data,threshold1,threshold2)
    mutate(data,posPSA=psa>=threshold1,posBP=bp>=threshold2) %>% group_by(advCancer,posPSA,posBP) %>% summarize(freq=n())
rTPF <- function(data) {
    a <- filter(data,advCancer & posPSA & posBP)$freq
    b <- filter(data,advCancer & !posPSA & posBP)$freq
    c <- filter(data,advCancer & posPSA & !posBP)$freq
    (a+b)/(a+c)
}
rFPF <- function(data) {
    e <- filter(data,!advCancer & posPSA & posBP)$freq
    f <- filter(data,!advCancer & !posPSA & posBP)$freq
    g <- filter(data,!advCancer & posPSA & !posBP)$freq
    (e+f)/(e+g)
}
rBiopsy <- function(data) 
    sum(filter(data,posBP)$freq)/sum(filter(data,posPSA)$freq)
RNGkind("Mersenne-Twister")
set.seed(12345)
d <- psaSimCor(10000,age=70,mubeta2.scale=2.1,rho=0.62)
plot(log(bp) ~ log(psa), data=d)
cor(subset(d,select=c(lpsa,lbp)))
## dAgg(d,3,3)
dAgg(d,3,3) %>% rTPF()
dAgg(d,3,3) %>% rFPF()

## uniroot1 <- uniroot(function(x) dAgg(d,3,x) %>% rBiopsy()-1, interval=c(1,20))
uniroot1 <- uniroot(function(x) dAgg(d,3,x) %>% rTPF()-1, interval=c(1,20))
uniroot1
## dAgg(d,3,uniroot1$root)
dAgg(d,3,uniroot1$root) %>% rTPF()
dAgg(d,3,uniroot1$root) %>% rFPF()
dAgg(d,3,uniroot1$root) %>% rBiopsy()

## Random draw from a bivariate normal distribution
rho <- 0.62
Sigma <- matrix(c(1,rho,rho,1),2)
A <- chol(Sigma)
z <- matrix(rnorm(2*1e5),nrow=1e5)
y <- cbind(z[,1],z[,1]*rho+z[,2]*sqrt(1-rho*rho))
## y <- z %*% A
cor(y)
apply(y,2,mean)
apply(y,2,sd)




lpsa1 <- psaSimCor(1e5,age=70,rho=0.0)
lpsa2 <- apply(lpsa1,1,mean)
## lpsa2 <- apply(lpsa1,1,function(x) sum(x*c(0.3,0.7)))
cor(lpsa1[,1],lpsa2) # cor>=0.71


plot(density(psaSim(1e4)$psa),xlim=c(0,20))
i <- 1
for (age in seq(55,80,by=5)) {
    lines(density(psaSim(1e5,age=age)$psa),col=i)
    i <- i+1
}





## The baseline FHCRC model assumes that PSA is an unbiased measure of the underlying diease process. The results here suggest that imprecision in the measure is less important than bias - and that PSA would need to be relatively biased to get the predicted change in biopsies from STHLM3.
## The main challenge now is that the FHCRC model was based on PCPT trial data which will not be available - nor, probably, will the bias be estimable from observed data.

## correlated betas
rho <- 0.5 # correlation
set.seed(12345+1)
beta0 <- correlatedValue(temp2$beta0,p$mubeta0,p$sebeta0, rho)
beta1 <- correlatedValue(temp2$beta1,p$mubeta1,p$sebeta1, rho)
beta2 <- pmax(0,correlatedValue(temp2$beta2,p$mubeta2[temp2$grade],p$sebeta2[temp2$grade], rho)) # should be a conditional distribution
biomarker2 <- exp(beta0+beta1*(temp2$age-35)+beta2*pos(temp2$age-35-temp2$t0)+rnorm(n,0,sqrt(p$tau2)))

## completely independent biomarker with more measurement error
set.seed(12345+1)
beta0 <- rnorm(n,p$mubeta0,p$sebeta0)
beta1 <- rnorm(n,p$mubeta1,p$sebeta1)
beta2 <- rnormPos(n,p$mubeta2[temp2$grade],p$sebeta2[temp2$grade])
biomarker2 <- exp(beta0+beta1*(temp2$age-35)+beta2*pos(temp2$age-35-temp2$t0)+rnorm(n,0,2*sqrt(p$tau2)))
plot(temp2$psa,biomarker2,log="xy")

set.seed(12345+2)
temp2 <- transform(temp2,
                   BBP=Z)
## STHLM3 simulation report
with(list(threshold=3.11),with(transform(temp2,BBPpos=(psa>=1 & BBP>=threshold),PSApos=(psa>=3)),
     cat(sprintf("
PSA+ & advanced:\t\t%i
PSA+ & cancer:\t\t\t%i
BBP+ & adv:\t\t\t%i
BBP+ & can:\t\t\t%i
BBP+ & PSA+:\t\t\t%i
BBP+ & PSA-:\t\t\t%i
BBP- & PSA+:\t\t\t%i
PSA+ | BBP+:\t\t\t%i
Prop reduction in biospies:\t%5.3f\n",
                 sum(PSApos & advanced),
                 sum(PSApos & cancer),
                 sum(BBPpos & advanced),
                 sum(BBPpos & cancer),
                 sum(BBPpos & PSApos),
                 sum(BBPpos & !PSApos),
                 sum(!BBPpos & PSApos),
                 sum(BBPpos | PSApos),
                 (sum(!BBPpos & PSApos) - sum(BBPpos & !PSApos))/
                       sum(BBPpos | PSApos))))))

logZ=rnorm(100000,0,0.1)
logpsa=logZ+rnorm(100000,0,sqrt(p$tau2))
Z=exp(logZ)
psa=exp(logpsa)
plot(density(logZ))
lines(density(logpsa),lty=2)
mean(logZ)
mean(logpsa)
mean(Z)
mean(psa)



## simulating sequentially from a bivariate normal distribution
set.seed(12345)
n <- 1e5
rho <- 0.6
sigma <- 2
u1 <- rnorm(n)
u2 <- rnorm(n)
x1 <- u1*sigma
x2 <- rho*u1*sigma+sqrt(1-rho^2)*u2*sigma
cor(x1,x2)
var(x1)
var(x2)

    
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
levelplot(pred~age*year,dx,subset=(age>=30), col.regions = rev(brewer.div(100)), aspect = "iso",
  xlab="Age (years)", ylab="Calendar period", ylim=c(1980,2050))
dev.off()
with(list(res=600),
       jpeg(file="~/work/levelplot-pc-incidence.jpg",height=5*res,width=5*res,res=res,
            quality=100))
levelplot(pred~age*year,dx,subset=(age>=30), col.regions = rev(brewer.div(100)), aspect = "iso",
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

require(microsimulation)
require(dplyr)
require(splines)
base <- callFhcrc(1e6,screen="noScreening",mc.cores=3,cohort=1970)
new <- callFhcrc(1e6,screen="twoYearlyScreen50to70",mc.cores=3,cohort=1970)
baserate <- predict(base,"cancerdeath") %>% filter(age>=50)
newrate <- predict(new,"cancerdeath") %>% filter(age>=50)
merged <- rbind(transform(baserate,group=0), transform(newrate,group=1)) %>% transform(age50=age-50)
fit1 <- glm(n ~ ns(age,5)+group:ns(age,5)+offset(log(pt)), data=merged, family=poisson)
RR <- predict(fit1,newdata=transform(baserate,group=1,pt=1,age50=age-50),type="response") /
predict(fit1,newdata=transform(baserate,group=0,pt=1,age50=age-50),type="response")

pdf("~/work/mortality_rate_reduction_twoYearly50to69.pdf",width=7,height=4)
par(mfrow=1:2)
with(predict(new,"incidence") %>% filter(age>=40),
     plot(age,rate*1e5,type="l",xlab="Age (years)",ylab="Prostate cancer incidence rate per 100,000",main="(a)"))
legend("topleft",legend=c("No screening","Screening"),lty=2:1, bty="n")
with(predict(base,"incidence") %>% filter(age>=40),
     lines(age,rate*1e5,lty=2))
plot(baserate$age, RR, type="l",ylab="Prostate cancer mortality rate ratio",xlab="Age (years)",
ylim=c(0.5,1),main="(b)")
dev.off()

with(predict(new,"cancerdeath") %>% filter(age>=40),
     plot(age,rate*1e5,type="l",xlab="Age (years)",ylab="Rate per 100,000"))
with(predict(base,"cancerdeath") %>% filter(age>=40),
     lines(age,rate*1e5,lty=2))



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
