library(igraphdata)
library(igraph)
data("karate")

k_edges <- as.data.frame(as_edgelist(karate))

write.csv(k_edges,
          "Data/karate_edge_list.csv",
          quote = FALSE,
          row.names = FALSE)
