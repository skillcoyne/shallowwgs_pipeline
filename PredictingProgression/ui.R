#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Barrett's Progression Prediction"),
  sidebarLayout(
    sidebarPanel(
      tags$h3('Patient Information'),
      numericInput("ageY", "Age at diagnosis:", min=1, max=100, value=60,step=1),
      sliderInput('beLength', 'Length of BE segment', min=1,max=30,round=1,value=4),
      radioButtons("sex", label = "Sex", choices = list("Male" = 'M', "Female" = 'F'), selected = 'M'),
      
      
      tags$h3('Patient Data'),
      fileInput("pathfile", "Optional: upload a file with per-sample pathology and/or p53 IHC (see example)",
                accept = c(
                  "text/tsv",
                  "text/csv",
                  "text/tab-separated-values,comma-separated-values,text/plain",
                  ".csv", ".tsv", ".txt")
      ),
      
      fileInput("file", "Choose Caveman copy number or segmented copy number file:",
                accept = c(
                  "text/tsv",
                  "text/csv",
                  "text/tab-separated-values,comma-separated-values,text/plain",
                  ".csv", ".tsv", ".txt")
      )
      # Input: Select separator ----
      #checkboxInput("header", "Header", TRUE),
      #checkboxInput("norm", "Model provided normalization (default: yes)", TRUE),
      
      #actionButton("action", label = "Submit")

    ),
  
    mainPanel(
      helpText(tags$b("Recommendations are evaluated with the following criteria:"), tags$br(),
               tags$ol(tags$li("Immediate RFA: HGD diagnosis or more than one consecutive high risk predictions."),
                       tags$li("Recheck 6-12 months: One high risk prediction or an aberrant p53 IHC"),
                       tags$li("Recheck endoscopy 12-24 months: One or more moderate risk predictions"),
                       tags$li("Regular surveillance 3-5 years: Two or more consecutive low risk predictions"))),
      
      tabsetPanel(
        type='tabs',
        tabPanel('Recommendations', plotOutput("plot"), plotOutput("riskTable")),
        tabPanel('Example data file',tableOutput("example")),
        tabPanel('Example pathology & p53 data file', tableOutput("examplep53"))
      )


      #uiOutput('samples')
      #plotOutput("plot"),
      #plotOutput("riskTable")
      
    )
  )
  
  
))
