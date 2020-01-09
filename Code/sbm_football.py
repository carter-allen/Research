# Import graph-tool
from graph_tool.all import *
from pylab import *

# Load football data
g = graph_tool.collection.data["football"]

# Fit degree corrected stochastic block model to football data
state = minimize_blockmodel_dl(g)

# Plot football network colored by group membership
state.draw(pos=g.vp.pos, output="Figures/football-sbm-fit.svg")

# Retrieve the fitted group membership
b = state.get_blocks()
# Print the membership of vertex 10
print(b[10])

# Access the matrix of edge counts between groups
e = state.get_matrix()
matshow(e.todense()) #matshow is from pylab
savefig("Figures/football-edge-counts.svg")


# Nested SBM
g = graph_tool.collection.data["celegansneural"]
state = minimize_nested_blockmodel_dl(g)
state.draw(output="Figures/celegans-hsbm-fit.svg")
state.print_summary()


# SBM with covariates
g = graph_tool.collection.konect_data["moreno_train"]
state = minimize_nested_blockmodel_dl(g,
                                      state_args=dict(recs=[g.ep.weight],
                                                      rec_types=["discrete-binomial"]))
state.draw(edge_color=g.ep.weight,
           ecmap=(matplotlib.cm.inferno, .6),
           eorder=g.ep.weight,
           edge_pen_width=prop_to_size(g.ep.weight, 1, 4, power=1),
           edge_gradient=[],
           output="Figures/moreno-train-wsbm.svg")
