## @knitr  ReadingAndManipulating
## ===========================================
## Setup
## ===========================================
rm(list=ls())
require(reshape)
require(ggplot2)
require(microsimulation)
require(data.table)

## ===========================================
## Helper functions
## ===========================================
eventRates <- function(obj,pattern="toClinicalDiagnosis",cumAge=NULL) {
    events <- data.table(obj$summary$events,key="event")
    pt <- data.table(obj$summary$pt,key="age")
    if (is.numeric(cumAge)){
        events[,age:=ifelse(age>cumAge,cumAge,age)]
        pt[,age:=ifelse(age>cumAge,cumAge,age)]
    }
    out <- with(merge(pt[,sum(pt),by=list(age,year)],
                      events[J(pattern),sum(n),by=list(age,year)], by=c("age","year")),
                transform(data.table(age=age, year=year,pt=V1.x,n=ifelse(is.na(V1.y),0.0,V1.y))[-1,],rate=n/pt))
    return(out)
}

## @knitr RunSim
## ===========================================
## Run simulation
## ===========================================    
n=1e6
noScreening <- callFhcrc(n=n, screen="noScreening")#, mc.cores=n.cores)
screenUptake <- callFhcrc(n=n, screen="screenUptake")#, mc.cores=n.cores)


## @knitr Preprocessing
## ===========================================
## Preprocessing
## ===========================================    

## Preparing the simulated data
simAllMortNoScreening <- data.table(cbind(eventRates(noScreening, pattern=c("toCancerDeath","toOtherDeath"), cumAge=100), origin="simNoScreening"))
simAllMortScreenUptake <- data.table(cbind(eventRates(screenUptake, pattern=c("toCancerDeath","toOtherDeath"), cumAge=100), origin="simScreenUptake"))
simPcMortNoScreening <- data.table(cbind(eventRates(noScreening, pattern="toCancerDeath", cumAge=85), origin="simNoScreening"))
simPcMortScreenUptake <- data.table(cbind(eventRates(screenUptake, pattern="toCancerDeath", cumAge=85), origin="simScreenUptake"))
simPcIncNoScreening <- data.table(cbind(eventRates(noScreening, pattern="toClinicalDiagnosis", cumAge=85), origin="simNoScreening"))
simPcIncScreenUptake <- data.table(cbind(eventRates(screenUptake, pattern=c("toClinicalDiagnosis","toScreenDiagnosis"), cumAge=85), origin="simScreenUptake"))

allCauseMortDt <- rbindlist(list(simAllMortNoScreening[,list(age,year,origin,rate)], 
                                 simAllMortScreenUptake[,list(age,year,origin,rate)]),
                            use.names=TRUE)

pcMortDt <- rbindlist(list(simPcMortNoScreening[,list(age,year,origin,rate)],
                           simPcMortScreenUptake[,list(age,year,origin,rate)]),
                           use.names=TRUE)

pcIncDt <- rbindlist(list(simPcIncNoScreening[,list(age,year,origin,rate)],  
                          simPcIncScreenUptake[,list(age,year,origin,rate)]),
                          use.names=TRUE)

## @knitr Plotting
## ===========================================
## All-cause mortality plotting
## ===========================================

#X11()
p <- ggplot(data = allCauseMortDt, aes(x=age, y=rate, colour=year))
p <- p + geom_line(aes(group = year)) + scale_y_log10()
p + facet_grid(.~ origin) + ggtitle("All-cause mortality, by data set")

#X11()
p <- ggplot(subset(allCauseMortDt,year>=1968 & year<=2013 & age >= 0), aes(age,rate))
p <- p + geom_line(aes(group=origin, colour=origin))
p + scale_y_log10() + facet_wrap( ~ year) + ggtitle("All-cause mortality, by all available calendar periods")

#X11()
p <- ggplot(subset(allCauseMortDt,year>=2005 & year<=2013 & age >= 0), aes(age,rate))
p <- p + geom_line(aes(group=origin, colour=origin))
p + scale_y_log10() + facet_wrap( ~ year) + ggtitle("All-cause mortality, by recent calendar periods")

## #X11()
## p <- ggplot(subset(allCauseMortDt,year>=1968 & year<=2013 & age >= 40), aes(year,rate))
## p <- p + geom_line(aes(group=origin, colour=origin))
## p + facet_wrap( ~ age) + ggtitle("All-cause mortality, by all available calendar periods")

#X11()
p <- ggplot(subset(allCauseMortDt,year>=1968 & year<=2013 & age >= 50 & age <= 70), aes(year,rate))
p <- p + geom_line(aes(group=origin, colour=origin))
p + facet_wrap( ~ age) + ggtitle("All-cause mortality, by ages 50-70")


## ===========================================
## Prostate cancer mortality
## ===========================================

#X11()
p <- ggplot(data = pcMortDt, aes(x=age, y=rate, colour=year))
p <- p+geom_line(aes(group = year)) + scale_y_log10()
p + facet_grid(.~ origin) + ggtitle("Prostate cancer mortality, by data set")

#X11()
p <- ggplot(subset(pcMortDt,year>=1997 & year<=2013 & age >= 30), aes(age,rate))
p <- p + geom_line(aes(group=origin, colour=origin)) + scale_y_log10()
p + facet_wrap( ~ year) + ggtitle("Prostate cancer mortality, by all available calendar periods")

#X11()
p <- ggplot(subset(pcMortDt,year>=2005 & year<=2013 & age >= 30), aes(age,rate))
p <- p + geom_line(aes(group=origin, colour=origin)) + scale_y_log10()
p + facet_wrap( ~ year) + ggtitle("Prostate cancer mortality, by recent calendar periods")

#X11()
ann_text <- data.frame(year=2005, rate=0.0025, lab = "Text", age = 85)
p <- ggplot(pcMortDt[year>=1997 & year<=2013 & is.element(age,c(seq(47,82,5),85))], aes(year,rate))
p <- p + geom_line(aes(group=origin, colour=origin))
p <- p + facet_wrap( ~ age) + ggtitle("Prostate cancer mortality, by all available ages")
p + geom_text(data = ann_text,label = "Cumulated at \nage 85+", size = 4)
#facet_wrap( ~ age, scales = "free")

#X11()
p <- ggplot(pcMortDt[year>=1997 & year<=2013 & is.element(age,seq(52,75,5))], aes(year,rate))
p <- p + geom_line(aes(group=origin, colour=origin))
p + facet_wrap( ~ age) + ggtitle("Prostate cancer mortality, by ages of interest")

## ===========================================
## Prostate cancer incidence
## ===========================================

#X11()
p <- ggplot(data = pcIncDt, aes(x=age, y=rate, colour=year))
p <- p + geom_line(aes(group = year)) + scale_y_log10()
p + facet_grid(.~ origin) + ggtitle("Prostate cancer incidence, by data set")

#X11()
p <- ggplot(subset(pcIncDt,year>=1970 & year<=2012 & age >= 30), aes(age,rate))
p <- p + geom_line(aes(group=origin, colour=origin)) + scale_y_log10()
p + facet_wrap( ~ year) + ggtitle("Prostate cancer incidence, by all available calendar periods")

#X11()
p <- ggplot(subset(pcIncDt,year>=2004 & year<=2012 & age >= 30), aes(age,rate))
p <- p + geom_line(aes(group=origin, colour=origin)) + scale_y_log10()
p + facet_wrap( ~ year) + ggtitle("Prostate cancer incidence, by recent calendar periods")

#X11()
ann_text <- data.frame(year=1980, rate=0.0025, lab = "Text", age = 85)
p <- ggplot(pcIncDt[year>=1970 & year<=2013 & is.element(age,c(seq(47,82,5),85))], aes(year,rate))
p <- p + geom_line(aes(group=origin, colour=origin))
p <- p + facet_wrap( ~ age) + ggtitle("Prostate cancer incidence, by all available ages")
p + geom_text(data = ann_text,label = "Cumulated at \nage 85+", size = 4)
#facet_wrap( ~ age, scales = "free")

#X11()
p <- ggplot(pcIncDt[year>=1970 & year<=2013 & is.element(age,seq(52,77,5))], aes(year,rate))
p <- p + geom_line(aes(group=origin, colour=origin))
p + facet_wrap( ~ age) + ggtitle("Prostate cancer incidence, by ages of interest")

#graphics.off()
