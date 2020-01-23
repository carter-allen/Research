# Python script to simulate data integration between two graphs
# Graph 1 (G1) is (large) literature mining graph
# Graph 2 (G2) is (small) gene expression graph
# Both graphs are sampled from an SBM with same probability parameter P
# Number of blocks B for G1 and G2 is fixed at B = 3
# Edges from G1 and G2 are combined to form a multigraph G

from graph_tool.all import *
import random as rand
import numpy as np
import pandas as pd
import os


n1 = 50 # number of nodes in G1
n2 = 50 # number of nodes in G2
N1 = list(range(0,n1)) # Node list for G1
N2 = list(rand.sample(N1,n2)) # Node list for G2
B = 3 # number of blocks in G1 and G2
b = list(range(0,B)) # list of possible groups (0,...,B-1)
P = np.array([[0.7000,0.0250,0.0250],
              [0.0250,0.7000,0.0250],
              [0.0250,0.0250,0.7000]]) # numpy probability matrix shared by G1 and G2
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
# Integrate G1 and G2
for e in G2.edges():
    G.add_edge(e.source(),e.target())
G.set_directed(False)

graph_draw(G,output = "Figures/integrated_graph.svg")

fit1 = minimize_blockmodel_dl(G1)
fit1.draw(output = "Figures/integrated_fit1.svg")

fit2 = minimize_blockmodel_dl(G2)
fit2.draw(output = "Figures/integrated_fit2.svg")

fit = minimize_blockmodel_dl(G)
fit.draw(output = "Figures/integrated_fit.svg")

# Simulation to compare integrated to non-integrated SBMs
I = 10 # number of simulations
n_G_correct = 0
n_G1_correct = 0
n_G2_correct = 0
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
    # Integrate G1 and G2
    for e in G2.edges():
        G.add_edge(e.source(),e.target())
    G.set_directed(False)
    
    fit1 = minimize_blockmodel_dl(G1,deg_corr = False)
    fit2 = minimize_blockmodel_dl(G2,deg_corr = False)
    fit = minimize_blockmodel_dl(G, deg_corr = False)
    bl1 = fit1.get_B()
    bl2 = fit2.get_B()
    bl = fit.get_B()
    print(bl,bl1,bl2)
    if bl1 == B: n_G1_correct = n_G1_correct + 1
    if bl2 == B: n_G2_correct = n_G2_correct + 1
    if bl == B: n_G_correct = n_G_correct + 1

# This approach of network integration seems not to work well.
