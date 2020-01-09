# Python script to simulate data integration between two graphs
# Graphs 1 and 2 (G1,G2) are two realizations of SBMs on same node set
# Both graphs are sampled from an SBM with same probability parameter P
# Number of blocks B for G1 and G2 is fixed at B = 3
# Edges from G1 and G2 are combined (edge union) to form a weighted simple graph G

from graph_tool.all import *
import random as rand
import numpy as np
from scipy import stats
import statistics
import pandas as pd
import os

n1 = 100 # number of nodes in G1
n2 = 100 # number of nodes in G2
N1 = list(range(0,n1)) # Node list for G1
N2 = N1 # Node list for G2
B = 3 # number of blocks in G1 and G2
b = list(range(0,B)) # list of possible groups (0,...,B-1)
P = np.array([[0.3000,0.1000,0.1000],
              [0.1000,0.3000,0.1000],
              [0.1000,0.1000,0.3000]]) # numpy probability matrix shared by G1 and G2
b1 = np.array(rand.choices(b, k = n1)) # sample block memberships for G1
b2 = b1 # block memberships for G2

I = 10 # number of simulations
p_fit_correct = np.empty(I)
p_fit1_correct = np.empty(I)
p_fit2_correct = np.empty(I)
n_G_correct = 0
n_G1_correct = 0
n_G2_correct = 0
sum_G = 0
sum_G1 = 0
sum_G2 = 0
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

    G = Graph(G1) # copy of G1
    for e in G2.edges():
        if G.edge(e.source(),e.target()) == None:
            G.add_edge(e.source(),e.target())

    fit = minimize_blockmodel_dl(G, deg_corr = False)
    fit1 = minimize_blockmodel_dl(G1,deg_corr = False)
    fit2 = minimize_blockmodel_dl(G2,deg_corr = False)

    bl1 = fit1.get_B()
    bl2 = fit2.get_B()
    bl = fit.get_B()

    sum_G = sum_G + bl
    sum_G1 = sum_G1 + bl1
    sum_G2 = sum_G2 + bl2

    print(i,bl,bl1,bl2)
    if bl1 == B: n_G1_correct = n_G1_correct + 1
    if bl2 == B: n_G2_correct = n_G2_correct + 1
    if bl == B: n_G_correct = n_G_correct + 1

    blks = np.empty(n1)
    blks1 = np.empty(n1)
    blks2 = np.empty(n2)
    for v in range(0,n1): # save node labels from each SBM to arrays
        blks[v] = fit.get_blocks()[v]
        blks1[v] = fit1.get_blocks()[v]
        blks2[v] = fit2.get_blocks()[v]

    n_fit_correct = np.empty(B)
    n_fit_correct1 = np.empty(B)
    n_fit_correct2 = np.empty(B)
    for a in range(0,B):
        n_fit_correct[a] = sum(blks[b1 == a] == stats.mode(blks[b1 == a])[0])
        n_fit_correct1[a] = sum(blks1[b1 == a] == stats.mode(blks1[b1 == a])[0])
        n_fit_correct2[a] = sum(blks2[b1 == a] == stats.mode(blks2[b1 == a])[0])
    p_fit_correct[i] = sum(n_fit_correct)/n1
    p_fit1_correct[i] = sum(n_fit_correct1)/n1
    p_fit2_correct[i] = sum(n_fit_correct1)/n1



