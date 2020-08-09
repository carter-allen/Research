plot <- function(x,...)
{
    UseMethod("plot",x)
}

plot.SBM <- function(fit)
{
    G = igraph::graph_from_adjacency_matrix(fit$A, mode = "undirected")
    plot(G, vertex.color = fit$z)
}