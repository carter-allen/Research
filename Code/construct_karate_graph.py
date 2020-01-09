import pandas as pd
from graph_tool.all import *

# Read karate edges as data frame
karate_edges = pd.read_csv("Data/karate_edge_list.csv")
V1 = karate_edges["V1"].tolist() # column 1
V2 = karate_edges["V2"].tolist() # column 2
actors = pd.unique(V1 + V2).tolist()  # unique set of actors (nodes)
# print(actors)
N = len(actors) # number of nodes in the network
ns = list(range(0,N))
k_dict = dict(zip(actors,ns))

g = Graph() # construct empty graph
g.add_vertex(N) # set nodes

# Populate graph with (directed) edges
# Note the karate network is built with directed edges
for i, r in karate_edges.iterrows(): # this is a pandas iterator for data frames
    rl = r.tolist()
    g.add_edge(g.vertex(k_dict[rl[0]]),g.vertex(k_dict[rl[1]]))

# Save simple network image
graph_draw(g,output="Figures/karate.svg")

# Modeling
fit1 = minimize_blockmodel_dl(g)
fit1.draw(output="Figures/karate_fit.svg")
fit_blks = fit1.get_blocks()
