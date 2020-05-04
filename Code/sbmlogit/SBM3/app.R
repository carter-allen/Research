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
setwd("~/Documents/School/Research/Code/sbmlogit/SBM3")

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
                        value = 2000,
                        step = 1000),
            actionButton("run",label = "Run Model"),
            actionButton("plot", label = "Plot Fit"),
            actionButton("show_posts", label = "Show Community Structure")
        ),
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("graphPlot"),
           plotOutput("fitPlot"),
           plotOutput("gammaPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    fit_model <- eventReactive(input$run,
                              {
                                  withProgress(message = "Fitting model", 
                                               value = 0,
                                               {
                                                   path = paste0("data/",input$dataset,".Rdata")
                                                   load(path)
                                                   fit = sbmlogit.mcmc(graph = g,
                                                                       alpha = input$K,
                                                                       nsamples = input$niter)
                                                   save(fit,file = "temp_fit.Rdata")
                                                   
                                               })
                                  
                              })
    
    show_fit <- eventReactive(input$plot,
                              {
                                  withProgress(message = "Plotting fit", 
                                               value = 0,
                                               {
                                                   load("temp_fit.Rdata")
                                                   f = plot_sbmlogit(fit) + ggtitle("Inferred partition:")
                                                   f
                                               })
                                  
                              })

    output$graphPlot <- renderPlot({
        path = paste0("data/",input$dataset,".Rdata")
        load(path)
        p = ggraph(g,layout = "kk") + 
            geom_node_point(size = 4) + 
            geom_edge_link(alpha = 0.50) + 
            theme_void() + 
            ggtitle("Observed graph:")
        p
    })
    
    output$fitPlot <- renderPlot({
        fit_model()
        show_fit()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
