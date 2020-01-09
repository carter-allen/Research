# Python script to simulate data integration between two graphs
# Graph 1 (G1) is (large) literature mining graph
# Graph 2 (G2) is (small) gene expression graph
# Both graphs are sampled from an SBM with same probability parameter P
# Number of blocks B for G1 and G2 is fixed at B = 3
# Edges from G1 and G2 are combined to form a weighted simple graph G

from graph_tool.all import *
import random as rand
import numpy as np
import pandas as pd
import os


n1 = 50 # number of nodes in G1
n2 = 20 # number of nodes in G2
N1 = list(range(0,n1)) # Node list for G1
N2 = list(rand.sample(N1,n2)) # Node list for G2
B = 3 # number of blocks in G1 and G2
b = list(range(0,B)) # list of possible groups (0,...,B-1)
P = np.array([[0.5000,0.0250,0.0250],
              [0.0250,0.5000,0.0250],
              [0.0250,0.0250,0.5000]]) # numpy probability matrix shared by G1 and G2
b1 = np.array(rand.choices(b, k = n1)) # sample block memberships for G1
b2 = np.array(rand.choices(b, k = n2)) # sample block memberships for G2

# Populate G1 with edges
G1 = Graph()
G1.add_vertex(n1)
for a in range(0,n1):
    for b in range(a,n1):
        ba = b1[a]
        bb = b1[b]
        p = P[ba,bb]
        if a == b: p = 0
        is_edge = np.random.binomial(1,p)
        if(is_edge): G1.add_edge(a,b)
G1.set_directed(False)

# Populate G2 with edges
G2 = Graph()
G2.add_vertex(n2)
for a in range(0,n2):
    for b in range(a,n2):
        ba = b2[a]
        bb = b2[b]
        p = P[ba,bb]
        if a == b: p = 0
        is_edge = np.random.binomial(1,p)
        if(is_edge): G2.add_edge(a,b)
G2.set_directed(False)

# Combine G1 and G2 into a multigraph
N = pd.unique(N1 + N2) # total node set for integrated graph (N should == N1)
n = len(N) # total number of nodes in G (n should == n1)
G = Graph(G1) # copy of G1
eprop = G.new_edge_property("int32_t")
G.edge_properties["weight"] = eprop
# Integrate G1 and G2
for e in G2.edges():
    if G.edge(e.source(),e.target()) == None:
        G.add_edge(e.source(),e.target())
        G.ep.weight[e] = 0
    else:
        G.ep.weight[e] = G.ep.weight[e] + 1
G.set_directed(False)

graph_draw(G,output = "Figures/integrated_graph.svg")

fit1 = minimize_blockmodel_dl(G1)
fit1.draw(output = "Figures/integrated_fit1.svg")

fit2 = minimize_blockmodel_dl(G2)
fit2.draw(output = "Figures/integrated_fit2.svg")

fit_w = minimize_blockmodel_dl(G,
                             state_args=dict(recs=[g.ep.weight],
                                             rec_types=["discrete-binomial"]))
fit_w.draw(output = "Figures/integrated_fit_w.svg")
fit = minimize_blockmodel_dl(G)
fit.draw(output = "Figures/integrated_fit_uw.svg")

# Simulation to compare integrated to non-integrated SBMs
I = 100 # number of simulations
n_G_correct = 0
n_G1_correct = 0
n_G2_correct = 0
n_Gw_correct = 0
n_Gm_correct = 0
sum_G = 0
sum_G1 = 0
sum_G2 = 0
sum_Gw = 0
sum_Gm = 0
for i in range(0,I):
    # Populate G1 with edges
    G1 = Graph()
    G1.add_vertex(n1)
    for a in range(0,n1):
        for b in range(a,n1):
            ba = b1[a]
            bb = b1[b]
            p = P[ba,bb]
            if a == b: p = 0
            is_edge = np.random.binomial(1,p)
            if(is_edge): G1.add_edge(a,b)
    G1.set_directed(False)

    # Populate G2 with edges
    G2 = Graph()
    G2.add_vertex(n2)
    for a in range(0,n2):
        for b in range(a,n2):
            ba = b2[a]
            bb = b2[b]
            p = P[ba,bb]
            if a == b: p = 0
            is_edge = np.random.binomial(1,p)
            if(is_edge): G2.add_edge(a,b)
    G2.set_directed(False)

    # Combine G1 and G2 into a multigraph
    N = pd.unique(N1 + N2) # total node set for integrated graph (N should == N1)
    n = len(N) # total number of nodes in G (n should == n1)
    G = Graph(G1) # copy of G1
    Gm = Graph(G1)
    eprop = G.new_edge_property("int32_t")
    G.edge_properties["weight"] = eprop
    # Integrate G1 and G2
    for e in G2.edges():
        if G.edge(e.source(),e.target()) == None:
            G.add_edge(e.source(),e.target())
            G.ep.weight[e] = 0
        else:
            G.ep.weight[e] = G.ep.weight[e] + 1
    Gm.add_edge(e.source(),e.target())
    G.set_directed(False)
    
    fit1 = minimize_blockmodel_dl(G1,deg_corr = False)
    fit2 = minimize_blockmodel_dl(G2,deg_corr = False)
    fit = minimize_blockmodel_dl(G, deg_corr = False)
    fit_w = minimize_blockmodel_dl(G,
                                   state_args=dict(recs=[g.ep.weight],
                                                   rec_types=["discrete-binomial"]))
    fit_m = minimize_blockmodel_dl(Gm, deg_corr = False)
    bl1 = fit1.get_B()
    bl2 = fit2.get_B()
    bl = fit.get_B()
    blw = fit_w.get_B()
    blm = fit_m.get_B()
    sum_G = sum_G + bl
    sum_G1 = sum_G1 + bl1
    sum_G2 = sum_G2 + bl2
    sum_Gm = sum_Gm + blm
    sum_Gw = sum_Gw + blw
    print(i,bl,bl1,bl2,blw,blm)
    if bl1 == B: n_G1_correct = n_G1_correct + 1
    if bl2 == B: n_G2_correct = n_G2_correct + 1
    if bl == B: n_G_correct = n_G_correct + 1
    if blw == B: n_Gw_correct = n_Gw_correct + 1
    if blm == B: n_Gm_correct = n_Gm_correct + 1



