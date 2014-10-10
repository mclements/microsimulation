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
## Read csv files
## ===========================================
myPath <- path.expand('~/src/ki/MortRateCalc')

## All-cause mortality plotting
popSwe_raw <- as.data.frame(read.csv2(file.path(myPath,"Alive_men_in_Sweden_by_age_and_year.csv"),
                                      header = TRUE, sep = "\t", dec = ",", fill = TRUE,
                                      comment.char = "", skip=2,fileEncoding="latin1"))

deathSwe_raw <- as.data.frame(read.csv2(file.path(myPath,"Dead_men_in_Sweden_by_age_and_year.csv"),
                                        header = TRUE, sep = ";", dec = ",", fill = TRUE,
                                        comment.char = "", skip=2,fileEncoding="latin1"))

popSthlm_raw <- as.data.frame(read.csv2(file.path(myPath,"Alive_men_in_Stockholm_by_age_and_year.csv"),
                                        header = TRUE, sep = "\t", dec = ",", fill = TRUE,
                                        comment.char = "", skip=2,fileEncoding="latin1"))

deathSthlm_raw <- as.data.frame(read.csv2(file.path(myPath,"Dead_men_in_Stockholm_by_age_and_year.csv"),
                                          header = TRUE, sep = "\t", dec = ",", fill = TRUE,
                                          comment.char = "", skip=2,fileEncoding="latin1"))

## Prostate cancer mortality
pcMortality_raw <- read.csv2(file.path(myPath,"ProstateCancerMortalitySocialstyrelsenStatistikdatabasen_2014-10-03 14_03_42.csv"), header = TRUE, sep = ";", dec = ".", comment.char = "", skip=1, na.strings=c("­­"))


## Prostate cancer incidence
pcIncidence_raw <- read.csv2(file.path(myPath,"ProstateCancerIncidenceSocialstyrelsenStatistikdatabasen_2014-10-03 13_56_35.csv"), header = TRUE, sep = ";", dec = ".", comment.char = "", skip=1, na.strings=c("­­"))


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

reshapeSocialStyrelsen <- function(data.in){
    tmp <- data.table(melt(data.in,id.vars = c("År", "Region")))
    data.out <- within(subset(tmp,variable!="Totalt"),
                       {origin <- Region
                        levels(origin) <- sub("Riket","SweStat",levels(origin))
                        levels(origin) <- sub("Stockholms län","SthlmStat",levels(origin))
                        rate <- value/100000 #data given in per 100 000
                        year <- År
                        ## age span mid-point imputation
                        age <- variable
                        age <- sub('X85.','X85.85', age)  #trick for the last value
                        age <- sub('X','', age)
                        age <- as.character(lapply(lapply(strsplit(age,"[.]"),as.numeric),mean))
                        rm(År,Region,value,variable)})
    data.out$age <- as.numeric(data.out$age)
    return(data.table(data.out))
    
}

reshapeSCB <- function(data.in){
    data.out <- within(subset(data.in,ålder!="Totalt"),{
        ålder <- sub(' år','',ålder)
        ålder <- sub('[+]','',ålder)
        age <- as.integer(ålder)
        origin <- region
        rm(ålder,region,kön)})
    data.out <- melt(data.out, id=c('age','origin'))
    data.out <- within(data.out,{
        year <- variable
        year <- sub('X','', year)
        levels(origin) <- sub("00 Riket","SweStat",levels(origin))
        levels(origin) <- sub("01 Stockholms län","SthlmStat",levels(origin))
        rm(variable)})
    data.out$year <- as.numeric(data.out$year)
    return(data.out)
}

rateSCB <- function(pop.in, mort.in){
    rate.out <- merge(pop.in,mort.in,by=c("age","year","origin"))
    rate.out <- within(rate.out, {
        rate=value.y/value.x
        rm(value.x, value.y)})
    return(data.table(rate.out))
}

## @knitr RunSim
## ===========================================
## Run simulation
## ===========================================    
n=1e7
#n.cores <- 3
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

## Preparing the SCB and socialstyrelsen stats
allMortSwe <- rateSCB(reshapeSCB(popSwe_raw),reshapeSCB(deathSwe_raw))
allMortSthlm <- rateSCB(reshapeSCB(popSthlm_raw),reshapeSCB(deathSthlm_raw))
pcMortality <- reshapeSocialStyrelsen(pcMortality_raw)
pcIncidence <- reshapeSocialStyrelsen(pcIncidence_raw)

allCauseMortDt <- rbindlist(list(allMortSthlm,
                                 allMortSwe,
                                 simAllMortNoScreening[,list(age,year,origin,rate)], 
                                 simAllMortScreenUptake[,list(age,year,origin,rate)]),
                            use.names=TRUE)

pcMortDt <- rbindlist(list(pcMortality,
                           simPcMortNoScreening[,list(age,year,origin,rate)],
                           simPcMortScreenUptake[,list(age,year,origin,rate)]),
                           use.names=TRUE)

pcIncDt <- rbindlist(list(pcIncidence,
                          simPcIncNoScreening[,list(age,year,origin,rate)],  
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
