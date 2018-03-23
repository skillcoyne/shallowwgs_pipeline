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
      ),
      
      
      tags$hr(),
      checkboxInput("header", "Header", TRUE),
      # Input: Select separator ----
      checkboxInput("norm", "Model provided normalization (default: yes)", TRUE),
      # Input: Select separator ----
      radioButtons("patients", "Samples are from -",
                   choices = c('Single patient' = "same",
                               'Multiple patients' = "mult"),
                   selected = "same")

    ),
    
    mainPanel(
      tableOutput("example"),
      tableOutput("examplep53"),
      #tableOutput("contents"),
      helpText(tags$b("Recommendations are evaluated with the following criteria:"), tags$br(),
      tags$ol(tags$li("Immediate RFA: HGD diagnosis or more than one consecutive high risk predictions."),
          tags$li("Recheck 6-12 months: One high risk prediction or an aberrant p53 IHC"),
          tags$li("Recheck endoscopy 12-24 months: One or more moderate risk predictions"),
          tags$li("Regular surveillance 3-5 years: Two or more consecutive low risk predictions"))),
      plotOutput("plot"),
      plotOutput("riskTable")
      
    )
  )
  
  
))
