#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
setwd("~/Documents/Research/Code/viz")
library(shiny)
library(igraph)
library(tidygraph)
library(ggraph)
library(scran)

load("D_exp_z.RData")
load("D_spa_z.RData")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Mouse Cortex scRNA-seq + Spatial Vizualization"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("w",
                        "Mixing weight",
                        min = 0,
                        max = 1,
                        value = 0)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        # x    <- faithful[, 2]
        # bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        # hist(x, breaks = bins, col = 'darkgray', border = 'white')
        
        w = input$w
        D = (1-w)*D1 + w*D2
        G <- buildKNNGraph(D, k = 31)
        
        g <- G %>%
            as_tbl_graph() 
        p <- ggraph(g,layout = "kk") + 
            geom_node_point()  
        p
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
