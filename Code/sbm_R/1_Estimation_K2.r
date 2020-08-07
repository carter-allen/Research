library(igraph) ; library(MCMCpack)
setwd("~/Documents/School/Research/Code/sbm_R")
rm(list = ls())

# input data
DATA = read.graph("karate.gml",format=c("gml"))
DATA

Adj_mat = as_adjacency_matrix(DATA)
Adj_mat ; class(Adj_mat)
A_mat = as.matrix(Adj_mat)

A_mat[1,1]
A_mat[1,2]

n_sample = dim(A_mat)[[1]]

Club = c(1,1,1,1,1, 
         1,1,1,1,2, 
         1,1,1,1,2, 
         2,1,1,2,1, 
         2,1,2,2,2, 
         2,2,2,2,2, 
         2,2,2,2) # 1: Mr. Hi's, 2:Officer's
length(Club) ; sum(Club==1)
Club_shape = rep(NA,n_sample)
Club_shape[which(Club==1)] = "circle" ; Club_shape[which(Club==2)] = "square"

# setting / starting values
K = 2
z_vec = sample(x=(1:K), size=n_sample, replace=T)
n_vec = rep(0,K)
for (k in 1:K){
	n_vec[k] = sum(z_vec==k)
}
pi_vec = n_vec / n_sample

P_mat = array(0.4,c(K,K))

# MCMC
n_iter = 1000 ; alpha_hyp = 0.5 ; beta1_hyp = beta2_hyp = 0.5 

draw_z = array(0,c(n_iter,n_sample))
draw_pi = array(0,c(n_iter,K))
draw_P = array(0,c(n_iter,K,K))
draw_logf_z_given_A = rep(0,n_iter)

for (i_iter in 1:n_iter){
	
	# Step 1. pi
	for (k in 1:K){
		n_vec[k] = sum(z_vec==k)
	}
	alpha_star = (n_vec + alpha_hyp) 
	pi_vec = rdirichlet(n=1, alpha=alpha_star)
	
	# Step 2. z
	for (i_sample in 1:n_sample){
		pi_star = rep(0,K)
		for (k in 1:K){
			pi_star[k] = pi_vec[k]
			for (j_sample in c(1:n_sample)[-i_sample]){
				pi_star[k] = pi_star[k] * (P_mat[k,z_vec[j_sample]])^A_mat[i_sample,j_sample] * (1.0-P_mat[k,z_vec[j_sample]])^(1-A_mat[i_sample,j_sample])
			} # for (j_sample in c(1:n_sample)[-i_sample])
		} # for (k in 1:K)
		pi_star = pi_star / sum(pi_star)
		z_vec[i_sample] = sample(x=(1:K), size=1, prob=pi_star)
	} # for (i_sample in 1:n_sample)
	
	# Step 3. P
	for (a in 1:K){
		for (b in a:K){
			sum_A_ij = sum_one_minus_A_ij = 0 
			for (i_sample in 1:(n_sample-1)){
				for (j_sample in (i_sample+1):n_sample){
					if ( (z_vec[i_sample]==a)&(z_vec[j_sample]==b) ){
						sum_A_ij = sum_A_ij + A_mat[i_sample,j_sample]
						sum_one_minus_A_ij = sum_one_minus_A_ij + (1-A_mat[i_sample,j_sample])
					} # if ( (z_vec[i_sample]==a)&(z_vec[j_sample]==b) )
				} # for (j_sample in (i_sample+1):n_sample)
			} # for (i_sample in 1:n_sample)
			beta1_star = beta1_hyp + sum_A_ij ; beta2_star = beta2_hyp + sum_one_minus_A_ij
			P_mat[a,b] = P_mat[b,a] = rbeta(n=1,beta1_star,beta2_star) 
		} # for (b in (a+1):K)
	} # for (a in 1:(K-1))

	# Calculating P(z|Data) to find MAP
	logf_z_given_A = 0
	for (i_sample in 1:n_sample){
		logf_z_given_A = logf_z_given_A + log(pi_vec[z_vec[i_sample]])
	} # for (i_sample in 1:n_sample)
	for (i_sample in 1:(n_sample-1)){
		for (j_sample in (i_sample+1):n_sample){
			logf_z_given_A = logf_z_given_A + A_mat[i_sample,j_sample]*log(P_mat[z_vec[i_sample],z_vec[j_sample]]) + (1.0-A_mat[i_sample,j_sample])*log(1.0-P_mat[z_vec[i_sample],z_vec[j_sample]])
		} # for (j_sample in (i_sample+1):n_sample)
	} # for (j_sample in (i_sample+1):n_sample)
	
	# Store values
	draw_z[i_iter,] = z_vec
	draw_pi[i_iter,] = pi_vec
	draw_P[i_iter,,] = P_mat
	draw_logf_z_given_A[i_iter] = logf_z_given_A
	
	# Print status
	if ( round(i_iter/100)==(i_iter/100) ){
		print( paste0("iter=",i_iter) )
		print( "pi_vec=" )
		print( round(pi_vec,3) )
	}
	
} # for (i_iter in 1:n_iter)

save.image("1_Estimation_K2.RData")

# Results
png(file="1_Estimation_K2a.png",width=800,height=500,pointsize=20)
par(mfrow=c(2,1),mai=c(0.8,0.8,0.4,0.4),family="serif",mgp = c(1.5, 0.5, 0)) # b l t r
plot(1:n_iter, draw_logf_z_given_A, type="l")
matplot(1:n_iter, draw_pi, type="l")
dev.off()

png(file="1_Estimation_K2b.png",width=1000,height=900,pointsize=20)
which_MAP = which.max(draw_logf_z_given_A)
MAP_grouping_result = draw_z[which_MAP,]
plot(DATA, vertex.color=MAP_grouping_result, vertex.shape=Club_shape )
dev.off()
