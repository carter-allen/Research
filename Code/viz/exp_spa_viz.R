# how to visualize expression and spatial data
setwd("~/Documents/Research/Code/viz")

# packages
library(igraph)
library(tidygraph)
library(ggraph)
library(scran)

# load sfp data (Giotto paper)
load("D_exp.RData")
load("D_spa.RData")
load("fit_A1.RData")
z1 = as.factor(fit_A1$z)

# normalize matrices
m1 <- mean(D1[lower.tri(D1)])
s1 <- sd(D1[lower.tri(D1)])
D1 <- (D1 - m1)/s1

m2 <- mean(D2[lower.tri(D2)])
s2 <- sd(D2[lower.tri(D2)])
D2 <- (D2 - m2)/s2

save(D1,file = "D_exp_z.RData")
save(D2,file = "D_spa_z.RData")

# integrate distance data
w = 0.2
D = (1-w)*D1 + w*D2
G <- buildKNNGraph(D, k = 31)

g <- G %>%
    as_tbl_graph() 
p <- ggraph(g,layout = "kk") + 
    geom_node_point()  
