#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(sbmlogit)
library(sbmlhelpers)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(igraphdata)
setwd("/home/carterallen/Documents/School/Research/Code/sbmlogit/SBM3")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("SBM Logit Interface"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("dataset",
                        label = "Select data:",
                        choices = list(
                            "Chesapeake (Lower)" = "ChesLower",
                            "Chesapeake (Middle)" = "ChesMiddle",
                            "Chesapeake (Upper)" = "ChesUpper",
                            "Karate" = "karate",
                            "Maspalomas" = "Maspalomas",
                            "U.K. Faculty" = "UKfaculty"
                        ),
                        selected = "karate"),
            tags$hr(),
            h4("Model parameters"),
            numericInput("K",
                        "Number of clusters:",
                        min = 1,
                        max = 5,
                        value = 2),
            sliderInput("niter",
                        "Number of MCMC iterations",
                        min = 1000,
                        max = 10000,
                        value = 2000),
            actionButton("run",label = "Run Model")
        ),
        # Show a plot of the generated distribution
        mainPanel(
           h3("Observed Data"),
           plotOutput("graphPlot"),
           h3("Inferred Partition"),
           plotOutput("fitPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    
    show_fit <- eventReactive(input$run,
                              {
                                  path = paste0("data/",input$dataset,".Rdata")
                                  load(path)
                                  fit = sbmlogit.mcmc(graph = g,
                                                      alpha = input$K,
                                                      nsamples = input$niter)
                                  f = plot_sbmlogit(fit)
                                  f
                              })

    output$graphPlot <- renderPlot({
        path = paste0("data/",input$dataset,".Rdata")
        load(path)
        p = ggraph(g,layout = "kk") + 
            geom_node_point(size = 4) + 
            geom_edge_link(alpha = 0.50) + 
            theme_void()
        p
    })
    
    output$fitPlot <- renderPlot({
        show_fit()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
