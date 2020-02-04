# function to plot SBM fit for sbmlogit package
library(ggraph)
library(tidygraph)
library(tidyverse)
source("~/Documents/School/Research/Code/sbmlogit/get_labels.R")

plot_sbmlogit <- function(fit, ground = NULL)
{
    # make sure graph is an igraph object
    if(class(fit$graph) != "igraph")
    {
        print("Error: provided graph is not an igraph object.")
        return(NULL) # return nothing
    }
    
    sigma = as.factor(get_labels(fit)) # returns the remapped posterior labels
    
    # if the user provided a ground truth variable, we want to map that to node shape
    # note that we assume the user provides a string containing the name of the 
    # ground truth cluster labeling in the as_tbl_graph() version of the igraph object
    if(!is.null(ground)) # if this parameter is null, then there is no ground truth
    {
        # if the ground parameter is not a charachter then just ignore it and
        # instead print the graph without any ground truth value mapped to the shapes
        if(!is.character(ground)) 
        {
            print("Error: ground truth parameter must be character. Ignoring ground truth parameter.")
            
            graph = fit$graph %>% # should be an igraph object
                as_tbl_graph() %>%
                activate(nodes) %>%
                mutate(sig = sigma,
                       eta = colMeans(fit$eta),
                       id = as.numeric(V(fit$graph)))
            
            p = ggraph(graph,layout = "kk") + 
                geom_edge_link(alpha = 0.5) +
                geom_node_point(aes(color = sig,
                                    size = eta)) + 
                geom_node_text(aes(label = id),
                               size = 3) + 
                theme_void()
        }
        # in the case the user provided a valid string parameter
        else
        {
            ground = fit$graph %>%
                as_tbl_graph() %>%
                activate(nodes) %>%
                pull(ground)
            
            graph = fit$graph %>% # should be an igraph object
                as_tbl_graph() %>%
                activate(nodes) %>%
                mutate(sig = sigma,
                       true = as.factor(ground),
                       eta = colMeans(fit$eta),
                       id = as.numeric(V(fit$graph)))
            
            p = ggraph(graph,layout = "kk") + 
                geom_edge_link(alpha = 0.5) +
                geom_node_point(aes(color = sig,
                                    shape = true,
                                    size = eta)) + 
                geom_node_text(aes(label = id),
                               size = 3) + 
                theme_void()
        }
    }
    # no ground truth parameter provided
    else 
    {
        graph = fit$graph %>% # should be an igraph object
            as_tbl_graph() %>%
            activate(nodes) %>%
            mutate(sig = sigma,
                   eta = colMeans(fit$eta),
                   id = as.numeric(V(fit$graph)))
        
        p = ggraph(graph,layout = "kk") + 
            geom_edge_link(alpha = 0.5) +
            geom_node_point(aes(color = sig,
                                size = eta)) + 
            geom_node_text(aes(label = id),
                           size = 3) + 
            theme_void()
    }
    
    return(p)
}