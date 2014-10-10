library(shiny)
library(BH)
library(Rcpp)
library(devtools)
#devtools::install_github('mclements/microsimulation@shinyApp')
#runGitHub("microsimulation", "mclements",subdir='shinyApp',ref="shinyApp") 
library(microsimulation)
library(parallel)
library(Epi)
library(sqldf)

# options(warn=-1)
# options(shiny.maxRequestSize=100*1024^2)


# Define server logic required to plot various variables
shinyServer(function(input, output) {
#   errind <- 0
  compliance <- 0.75
  participation <- 1.0

  ## ===============================================================
  ## population description
  plotLexis <- function(obj) {
    history <- obj$lifeHistories
    param <- obj$parameters
    tab <- sqldf("select t1.*, ageAtCancerDiagnosis, cohort, t0 from (select id, end as ageAtDeath, (event='toCancerDeath') as cancerDeath from history where event in ('toOtherDeath','toCancerDeath')) as t1  inner join param as p on p.id=t1.id left join (select id, end as ageAtCancerDiagnosis from history where event in ('toClinicalDiagnosis','toScreenDiagnosis')) as t2 on t1.id=t2.id")
    lexis1 <- Lexis(entry=list(coh=cohort,age=0),exit=list(coh=cohort+ageAtDeath,age=ageAtDeath),
                    data=tab)
    plot(lexis1, xlab="Calendar period", ylab="Age (year)", ylim=c(0,100), asp=1, main = 'Life Histories',
         cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    with(subset(tab,!is.na(ageAtCancerDiagnosis)),
         points(cohort+ageAtCancerDiagnosis,ageAtCancerDiagnosis,pch=19,cex=0.4,col="red"))
    with(subset(tab,t0+35<ageAtDeath),
         points(cohort+t0+35,t0+35,pch=19,cex=0.4,col="blue"))
    legend("topleft",legend=c("Latent cancer onset","Cancer diagnosis"),
           pch=19,col=c("blue","red"),bty="n",cex=1)
  }
  
  ## ===============================================================
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
  
  
  eventRates <- function(obj,pattern="Diagnosis") {
    stopifnot(require(sqldf))
    ev <- data.frame(event=grep(pattern,levels(obj$summary$events$event),value=TRUE))
    pt <- obj$summary$pt
    events <- obj$summary$events
    sqldf("select year, sum(pt) as pt, sum(n) as n, sum(rate*wt) as rate from (select cohort+age as year, age, pt, coalesce(n,0.0) as n, coalesce(n,0.0)/pt as rate from (select cohort, age, sum(pt) as pt from pt group by cohort, age) as t1 natural left outer join (select cohort, age, sum(n) as n from events natural join ev group by cohort, age) as t2) as main natural join w where year>=1990 and year<2030 group by year")
  }
  
  plotEvents <- function(obj,pattern, ylab="%", main=NULL) {
    with(eventRates(obj,pattern),
         plot(year, rate *100, type="l",ylim=c(0,max(eventRates(obj,pattern)$rate*100)),
              xlab="Calendar period", ylab=ylab, main=main,
              cex.main=1.5, cex.lab=1.5, cex.axis=1.5))
  }
  
  ## Data to be displayed at a valid simulation
  dataInput <- reactive({  
    if(input$get == 0) return(NULL)
    
    SimPop <<- isolate({  
      input$n
    })
    
    SceName <<- isolate({          
      getSceName <- function(prot){
        scenario_list <- list("No screening" = "noScreening", 
                              "Opportunistic screening" = "screenUptake", 
                              "Göteborg protocol (2+2)" = "stockholm3_goteborg",
                              "Risk-stratified protocol (4+8)" = "stockholm3_risk_stratified")
        names(scenario_list)[scenario_list %in% prot]}
      getSceName(prot=input$protocol)
    })
    
    IncludeLexis <<- isolate({
      input$IncludeLexis      
    })
    
    IncludeRates <<- isolate({TRUE})
    
    out <<- isolate({ 
      if (input$IncludeLexis){
        nLifeHistories=input$n
      }else{
        nLifeHistories=1
      }
      # Bug when using multiple cores?
      #mc.cores <- min(input$n, detectCores( logical = FALSE))
      callFhcrc(input$n, screen=input$protocol, nLifeHistories, mc.cores=1)      
    })
    
    InputText <<- isolate({paste('Results from simulating <span style="color:blue">', SceName,' </span> with a population of <span style="color:blue">', SimPop , ' individuals</span>.')})
  })
  
  
  ## Data to be displayed at an input error
  dataError <- reactive({  
    if(input$get == 0) return(NULL)
    
    IncludeLexis <<- isolate({FALSE})
    
    IncludeRates <<- isolate({FALSE})
    
    out <<- isolate({NULL})
    
    SimPop <<- isolate({input$n})
    
    SceName <<- isolate({NULL})
    
    InputText <<- isolate({paste('<span style="color:red">', SimPop,' </span> is not allowed as a population size.')})
  })
    
    ## Condition to run sim and update output
    finalInput <- reactive({
        ifelse(isolate({input$n<=100000 & input$n>0 & input$n%%1==0}), return(dataInput()), return(dataError()))
    })
  
  
  output$CurrentSimulation <- renderText({
    if(input$get == 0) return(NULL)
    finalInput()
    InputText
  })

  output$PlotsOut <- renderPlot({
    if(input$get == 0) return(NULL)
    
    
    if (IncludeRates){
    par(mfrow=c(2,2))
    rate_data <- out
    plotEvents(rate_data,"^toScreen$", main = 'PSA uptake')
    plotEvents(rate_data,"Biopsy", main = 'Biopsy uptake')
    plotEvents(rate_data,"Diagnosis", main="Prostate cancer incidence")
    plotEvents(rate_data,"^toCancerDeath$", main="Prostate cancer mortality")
    }
  })
  
  
  output$lexis <- renderPlot({
    if(input$get == 0) return(NULL)
    
    if (IncludeLexis){
      plotLexis(out)
    }
  })
})
