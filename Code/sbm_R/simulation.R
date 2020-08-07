library(MCMCpack)
library(sbmlhelpers)
library(sbmlogit)
library(tidygraph)
library(ggraph)
# setwd("~/Documents/School/Research/Code/sbm_R")
rm(list = ls())
#set.seed(1801)
Rcpp::sourceCpp("sample_SBM_fast.cpp")
Rcpp::sourceCpp("update_z_2.cpp")
Rcpp::sourceCpp("update_P.cpp")
Rcpp::sourceCpp("logf_z_given_A.cpp")

# GENERATE DATA
# Symmetric, non-reflexive 
n <- 500 # number of nodes
K <- 3
pi <- rep(1/K,K)
z <- sample(1:K,
            size = n,
            replace = TRUE,
            prob = pi)
p_in <- 0.50
p_out <- 0.05
P <- matrix(p_out,
            nrow = K,
            ncol = K)
diag(P) <- p_in
# A <- matrix(0, nrow = n, ncol = n)
# 
# for(i in 1:n)
# {
#     for(j in 1:i)
#     {
#         A[i,j] <- A[j,i] <- rbinom(1,1,P[z[i],z[j]]) # sample edges according to P and z
#     }
# }
# diag(A) <- 0
A <- sample_SBM_fast(z,P)
#G <- graph_from_adjacency_matrix(A,mode = "undirected")
#plot(G)

# Set priors
a0 = n/K # Dirichlet prior for pi
b10 = b20 = 0.5 # beta (uniform) priors for P elements

# Initialize parameters
K0 = 3 # putative number of clusters
zs = sample(1:K0, 
            size=n, 
            replace=TRUE) # initial cluster allocations
ns = table(zs) # initial cluster counts
pis = ns/n # initial cluster proportions
Ps = array(0.4,c(K0,K0)) # initial connectivity matrix

# MCMC settings and storage
n_iter = 1000 # number of iterations
burn = 100 # burn in 
n_sim = n_iter - burn # number of stored simulations
Z = array(0,c(n_sim,n)) # storage for cluster assignments
PI = array(0,c(n_sim,K0)) # storage for cluster proportions
PM = array(0,c(n_sim,K0,K0)) # storage for community connection params
draw_logf_z_given_A = rep(0,n_sim) # P(z|Data) for MAP estimate of z

start.time<-proc.time()
for (i in 1:n_iter){
    
    # Step 1. pi
    ns = table(zs)
    as = ns + a0
    pis = rdirichlet(n=1, as)
    print(ns)
    
    # Step 2. z
    # for (i_sample in 1:n){ # loop through all nodes
    #     pi_star = rep(0,K0)
    #     for (k in 1:K0){ # loop through all cluster
    #         pi_star[k] = pis[k]
    #         for (j_sample in c(1:n)[-i_sample]){ # loop through all nodes except i_sample
    #             pi_star[k] = pi_star[k] * 
    #                 (Ps[k,zs[j_sample]])^A[i_sample,j_sample] * 
    #                 (1.0-Ps[k,zs[j_sample]])^(1-A[i_sample,j_sample]) 
    #         } # for (j_sample in c(1:n_sample)[-i_sample])
    #     } # for (k in 1:K)
    #     pi_star = pi_star / sum(pi_star)
    #     # print(pi_star)
    #     zs[i_sample] = sample(x=(1:K0), size=1, prob=pi_star)
    # } # for (i_sample in 1:n_sample)
    zs = update_z(zs,A,Ps,pis,1:K0)
    
    # Step 3. P
    # for (a in 1:K0){
    #     for (b in a:K0){
    #         sum_A_ij = sum_one_minus_A_ij = 0
    #         for (i_sample in 1:(n-1)){
    #             for (j_sample in (i_sample+1):n){
    #                 if ( (zs[i_sample]==a)&(zs[j_sample]==b) ){
    #                     sum_A_ij = sum_A_ij + A[i_sample,j_sample]
    #                     sum_one_minus_A_ij = sum_one_minus_A_ij + (1-A[i_sample,j_sample])
    #                 } # if ( (z_vec[i_sample]==a)&(z_vec[j_sample]==b) )
    #             } # for (j_sample in (i_sample+1):n_sample)
    #         } # for (i_sample in 1:n_sample)
    #         beta1_star = b10 + sum_A_ij ; beta2_star = b20 + sum_one_minus_A_ij
    #         #print(paste(a,b,sum_A_ij,sum_one_minus_A_ij))
    #         Ps[a,b] = Ps[b,a] = rbeta(n=1,beta1_star,beta2_star)
    #     } # for (b in (a+1):K)
    # } # for (a in 1:(K-1))
    Ps = update_P(A,zs,K0,b10,b20)
    
    # Calculating P(z|Data) to find MAP
    # logf_z_given_A = 0
    # for (i_sample in 1:n){
    #     logf_z_given_A = logf_z_given_A + log(pis[zs[i_sample]])
    # } # for (i_sample in 1:n_sample)
    # for (i_sample in 1:(n-1)){
    #     for (j_sample in (i_sample+1):n){
    #         logf_z_given_A = logf_z_given_A + 
    #             A[i_sample,j_sample]*log(Ps[zs[i_sample],zs[j_sample]]) + 
    #             (1.0-A[i_sample,j_sample])*log(1.0-Ps[zs[i_sample],zs[j_sample]])
    #     } # for (j_sample in (i_sample+1):n_sample)
    # } # for (j_sample in (i_sample+1):n_sample)
    lz = logf_z_given_A(zs,A,pis,Ps)
    
    # Store values
    if(i > burn)
    {
        j = i - burn
        Z[j,] = zs
        PI[j,] = pis
        PM[j,,] = Ps
        draw_logf_z_given_A[j] = lz
    }
    
    # Print status
    print(paste("Iteration",i))
    #print(pis)
    if ( round(i/100)==(i/100) ){
        print( paste0("iter=",i) )
        print( "pi_vec=" )
        print( round(pis,3) )
    }
    
}
run.time<-proc.time()-start.time
print(paste("Finished MCMC after",run.time[1],"seconds"))

# Results
png(file="simulationa.png",width=800,height=500,pointsize=20)
par(mfrow=c(2,1),mai=c(0.8,0.8,0.4,0.4),family="serif",mgp = c(1.5, 0.5, 0)) 
plot(1:n_sim, draw_logf_z_given_A, type="l")
matplot(1:n_sim, PI, type="l")
dev.off()

png(file="simulationb.png",width=1000,height=900,pointsize=20)
which_MAP = which.max(draw_logf_z_given_A)
MAP_grouping_result = Z[which_MAP,]
#plot(G, vertex.color=MAP_grouping_result )
dev.off()

# fit sbmlogit for comparison
#fit_sbmlogit <- sbmlogit.mcmc(G,alpha = K0,nsamples = n_sim)
#z_sbmlogit <- get_labels(fit_sbmlogit)
#plot(G, vertex.color=z_sbmlogit )
#table(z_sbmlogit)
table(MAP_grouping_result)

#fit_MCMC <- fit_SBM_MCMC(A,K = 4)
table(MAP_grouping_result,zs)
