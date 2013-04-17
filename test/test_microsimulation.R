## try(detach("package:microsimulation", unload=TRUE))
## require(microsimulation)
## microsimulation:::.testPackage()

options(width=110)
require(microsimulation)
n <- 1e7
noScreening <- callFhcrc(n,screen="noScreening")
## "screenUptake", "stockholm3_goteborg", "stockholm3_risk_stratified"
uptake <- callFhcrc(n,screen="screenUptake")
goteborg <- callFhcrc(n,screen="stockholm3_goteborg",)
riskStrat <- callFhcrc(n,screen="stockholm3_risk_stratified")

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
  ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
  pt <- obj$summary$pt
  events <- obj$summary$events
  sqldf("select year, sum(pt) as pt, sum(n) as n, sum(rate*wt) as rate from (select cohort+age as year, age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select cohort, age, sum(pt) as pt from pt group by cohort, age) as t1 natural left outer join (select cohort, age, sum(n) as n from events natural join ev group by cohort, age) as t2) as main natural join w where year>=1990 and year<2020 group by year")
}

with(eventRates(noScreening),
     plot(year, rate, type="l",ylim=c(0,0.01),
          xlab="Age (years)", ylab="Rate"))
with(eventRates(uptake), lines(year, rate, col="red"))
with(eventRates(goteborg), lines(year, rate, col="green"))
with(eventRates(riskStrat), lines(year, rate, col="blue"))

## Plot of the cohorts over the Lexis diagram
plot(c(1900,2020),c(0,100),type="n",xlab="Calendar period",ylab="Age (years)")
polygon(c(1900,1970,1970+50,1970+50,1970+50-20,1900),
        c(0,0,50,100,100,0))
polygon(c(1990,2020,2020,1990),
        c(50,50,80,80),
        lty=2)



pdf("~/work/screening-comparison-20130222.pdf")
layout(matrix(1:4,nrow=2,byrow=TRUE))
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None verus two-yearly screening"))
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


### FHCRC model ###
options(width=110)
require(microsimulation)
n <- 1e5
temp <- callFhcrc(n,screen="noScreening")
temp2 <- callFhcrc(n,screen="twoYearlyScreen50to70")
temp4 <- callFhcrc(n,screen="fourYearlyScreen50to70")
temp50 <- callFhcrc(n,screen="screen50")
temp60 <- callFhcrc(n,screen="screen60")
temp70 <- callFhcrc(n,screen="screen70")

## rate calculations
require(sqldf)
eventRates <- function(obj,pattern="Diagnosis") {
  ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
  pt <- obj$summary$pt
  events <- obj$summary$events
  sqldf("select age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select age, sum(pt) as pt from pt group by age) as t1 natural left outer join (select age, sum(n) as n from events natural join ev group by age) as t2")
}

pdf("~/work/screening-comparison-20130222.pdf")
layout(matrix(1:4,nrow=2,byrow=TRUE))
with(eventRates(temp), plot(age, rate, type="l",
                            xlab="Age (years)", ylab="Rate", main="None verus two-yearly screening"))
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
