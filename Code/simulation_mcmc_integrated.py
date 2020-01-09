# Python script to simulate data integration between two graphs
# Graphs 1 and 2 (G1,G2) are two realizations of SBMs on same node set
# Both graphs are sampled from an SBM with same probability parameter P
# Number of blocks B for G1 and G2 is fixed at B = 3
# Edges from G1 and G2 are combined (edge union) to form a weighted simple graph G
# MCMC is used to get posterior marginals of edge marginals

from graph_tool.all import *
import random as rand
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import statistics
import pandas as pd
import os
my_colors = ['#BBBAD9','#91D9CC','#F27D72']

n1 = 100 # number of nodes in G1
n2 = 100 # number of nodes in G2
N1 = list(range(0,n1)) # Node list for G1
N2 = N1 # Node list for G2
B = 3 # number of blocks in G1 and G2
b = list(range(0,B)) # list of possible groups (0,...,B-1)
P = np.array([[0.50,0.10,0.10],
              [0.10,0.50,0.10],
              [0.10,0.10,0.50]]) # numpy probability matrix shared by G1 and G2
b1 = np.array(rand.choices(b, k = n1)) # sample block memberships for G1
b2 = b1 # block memberships for G2
blk_colors1 = ['empty'] * n1
blk_colors2 = ['empty'] * n2
for i in range(n1):
    blk_colors1[i] = my_colors[b1[i]]
    blk_colors2[i] = my_colors[b2[i]]

# Intialize G1 and set vertices
G1 = Graph() # Init
G1.add_vertex(n1) # Set vertices
G1_labels = G1.new_vertex_property("int32_t") # New vertex property for labels
G1_colors = G1.new_vertex_property("string") # New vertex property for colors
G1.vp.labels = G1_labels # Internalize vertex property
G1.vp.colors = G1_colors # Internalize color property
for v1 in range(n1): # Assign labels according to b1
    G1.vp.labels[v1] = b1[v1]
    G1.vp.colors[v1] = blk_colors1[v1]
# Populate G1 with edges according to an SBM(n1,P,b1)
for a in range(0,n1):
    for b in range(a,n1):
        ba = b1[a]
        bb = b1[b]
        p = P[ba,bb]
        if a == b: p = 0
        is_edge = np.random.binomial(1,p)
        if(is_edge): G1.add_edge(a,b)
G1.set_directed(False)
graph_draw(G1,
           vertex_text = G1.vp.labels,
           vertex_fill_color = G1.vp.colors,
           edge_color = "#d8dee8",
           output = "Figures/G1.svg") # Save observed graph viz

# Intialize G2 and set vertices
G2 = Graph() # Init
G2.add_vertex(n2) # Set vertices
G2_labels = G2.new_vertex_property("int32_t") # New vertex property
G2_colors = G2.new_vertex_property("string") # New vertex property for colors
G2.vp.labels = G2_labels # Internalize vertex property
G2.vp.colors = G2_colors # Internalize color property
for v2 in range(n2): # Assign labels according to b1
    G2.vp.labels[v2] = b2[v2]
    G2.vp.colors[v2] = blk_colors2[v2]
# Populate G2 with edges according to an SBM(n2,P,b2)
for a in range(0,n2):
    for b in range(a,n2):
        ba = b2[a]
        bb = b2[b]
        p = P[ba,bb]
        if a == b: p = 0
        is_edge = np.random.binomial(1,p)
        if(is_edge): G2.add_edge(a,b)
G2.set_directed(False)
graph_draw(G2,
           vertex_text = G2.vp.labels,
           vertex_fill_color = G2.vp.colors,
           edge_color = "#d8dee8",
           output = "Figures/G2.svg") # Save observed graph viz

G = Graph(G1) # copy of G1
for e in G2.edges():
    if G.edge(e.source(),e.target()) == None:
        G.add_edge(e.source(),e.target())
graph_draw(G,
           vertex_text = G.vp.labels,
           vertex_fill_color = G.vp.colors,
           edge_color = "#d8dee8",
           output = "Figures/G.svg") # Save observed graph viz

# Initialize block states of each network
fit = minimize_blockmodel_dl(G, deg_corr = False)
fit1 = minimize_blockmodel_dl(G1,deg_corr = False)
fit2 = minimize_blockmodel_dl(G2,deg_corr = False)

# Equilibrate MCMC chains
# mcmc_equilibrate(fit,wait = 1000,mcmc_args = dict(niter = 10))
# mcmc_equilibrate(fit1,wait = 1000,mcmc_args = dict(niter = 10))
# mcmc_equilibrate(fit2,wait = 1000,mcmc_args = dict(niter = 10))

# Collect and plot marginals for each graph
pv = None
def collect_marginals(s):
    global pv
    b = perfect_prop_hash([s.b])[0]
    pv = s.collect_vertex_marginals(pv, b=b)

mcmc_equilibrate(fit,
                 force_niter=10000,
                 mcmc_args=dict(niter=10),
                 callback=collect_marginals)
fit.draw(vertex_shape="pie",
         vertex_pie_fractions=pv,
         vertex_pie_colors = my_colors,
         vertex_text = G.vp.labels,
         edge_gradient=None,
         edge_color = "#d8dee8",
         output="Figures/mcmc_fit.svg")

pv = None
mcmc_equilibrate(fit1,
                 force_niter=10000,
                 mcmc_args=dict(niter=10),
                 callback=collect_marginals)
fit1.draw(vertex_shape="pie",
          vertex_pie_fractions=pv,
          vertex_pie_colors = my_colors,
          vertex_text = G1.vp.labels,
          edge_gradient=None,
          edge_color = "#d8dee8",
          output="Figures/mcmc_fit1.svg")

pv = None
mcmc_equilibrate(fit2,
                 force_niter=10000,
                 mcmc_args=dict(niter=10),
                 callback=collect_marginals)
fit2.draw(vertex_shape="pie",
          vertex_pie_fractions=pv,
          vertex_pie_colors = my_colors,
          vertex_text = G2.vp.labels,
          edge_gradient=None,
          edge_color = "#d8dee8",
          output="Figures/mcmc_fit2.svg")

# Collect posterior probability of B for each graph
def collect_num_groups(s):
    B = s.get_nonempty_B()
    h[B] += 1
h = np.zeros(G.num_vertices() + 1)
# Now we collect the marginals for exactly 100,000 sweeps, at
# intervals of 10 sweeps:
mcmc_equilibrate(fit,
                 force_niter=10000,
                 mcmc_args=dict(niter=10),
                 callback=collect_num_groups)
h_fit = h
plt.clf()
plt.bar(list(range(B-2,B+2)),h_fit[list(range(B-2,B+2))])
plt.xticks(list(range(B-2,B+2)),list(range(B-2,B+2)))
plt.savefig("Figures/fit_mcmc_B.svg")

h = np.zeros(G1.num_vertices() + 1)
# Now we collect the marginals for exactly 100,000 sweeps, at
# intervals of 10 sweeps:
mcmc_equilibrate(fit1,
                 force_niter=10000,
                 mcmc_args=dict(niter=10),
                 callback=collect_num_groups)
h_fit1 = h
plt.clf()
plt.bar(list(range(B-2,B+2)),h_fit1[list(range(B-2,B+2))])
plt.xticks(list(range(B-2,B+2)),list(range(B-2,B+2)))
plt.savefig("Figures/fit1_mcmc_B.svg")

h = np.zeros(G2.num_vertices() + 1)
# Now we collect the marginals for exactly 100,000 sweeps, at
# intervals of 10 sweeps:
mcmc_equilibrate(fit2,
                 force_niter=10000,
                 mcmc_args=dict(niter=10),
                 callback=collect_num_groups)
h_fit2 = h
plt.clf()
plt.bar(list(range(B-2,B+2)),h_fit2[list(range(B-2,B+2))])
plt.xticks(list(range(B-2,B+2)),list(range(B-2,B+2)))
plt.savefig("Figures/fit2_mcmc_B.svg")
