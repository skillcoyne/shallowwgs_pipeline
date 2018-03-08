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
      checkboxInput("norm", "Model provided norm", TRUE)

    ),
    
    
    
    mainPanel(
      tableOutput("contents"),
      plotOutput("plot")
      
    )
  )
  
  
))
