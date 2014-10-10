# Run with
# runApp("/home/andkar/src/ki/microsimulation_testing/shinyapp")

library(shiny)

# options(warn=-1)

# Code to make a message that shiny is loading
# Make the loading bar
loadingBar <- tags$div(class="progress progress-striped active",
                       tags$div(class="bar", style="width: 100%;"))
# Code for loading message
loadingMsg <- tags$div(class="modal", tabindex="-1", role="dialog", 
                       "aria-labelledby"="myModalLabel", "aria-hidden"="true",
                       tags$div(class="modal-header",
                                tags$h3(id="myModalHeader", "Running simulation...")),
                       tags$div(class="modal-footer",
                                loadingBar))

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Microsimulation demonstration"),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    numericInput("n", "Population size (1-100000):",
                10, min=1, max=100000),
    selectInput("protocol", "Screening protocol:",
                list("No screening" = "noScreening", 
                     "Opportunistic screening" = "screenUptake", 
                     "Göteborg protocol (2+2)" = "stockholm3_goteborg",
                     "Risk-stratified protocol (4+8)" = "stockholm3_risk_stratified")),
   actionButton("get", "Start simulation"),
   checkboxInput("IncludeLexis", "Include lexis diagram (takes a little longer).", FALSE)    
  ),
  
  # Show the caption and plot of the requested variables
  mainPanel(
    h3(htmlOutput("CurrentSimulation")),
    
    conditionalPanel( condition="input.IncludeLexis!=0 " ,
                      wellPanel(
                        h4("Lexis diagram of the individual life histories"),
                        helpText("The lexis diagram is informative when simulating a low number of individuals. It shows the individual life histories, including cancer onset and diagnosis."),
                        plotOutput("lexis")
                      )
    ),
    wellPanel(
      h4("Population rates"),
      helpText("The population rates are informative when simulating a high number of individuals. They show the PSA uptake, biopsy uptake, prostate cancer incidence (diagnosis) and prostate cancer mortality (deaths) expressed in percent per person-year. The screening effect of the implemented screening protocols is in addition to the oppertunistic (background) screening and the periodic pattern is the re-screening interval."),
      plotOutput("PlotsOut")
    ),
      
      ### Loading bar
      conditionalPanel(paste("$('html').hasClass('shiny-busy')"),
                       loadingMsg)
  )
))


