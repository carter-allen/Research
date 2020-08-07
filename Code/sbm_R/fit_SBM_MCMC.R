# Use Gibbs sampling to fit SBM
# For undirected, non-reflective graphs
fit_SBM_MCMC <- function(A,
                         K,
                         a0 = 0.5,
                         b10 = 0.5,
                         b20 = 0.5,
                         n_iter = 1000,
                         burn = 100)
{
    # Initialize parameters
    K0 = K # putative number of clusters
    zs = sample(1:K0, 
                size=n, 
                replace=TRUE) # initial cluster allocations
    ns = table(zs) # initial cluster counts
    pis = ns/n # initial cluster proportions
    Ps = array(0.4,c(K0,K0)) # initial connectivity matrix
    
    # MCMC settings and storage
    n_iter = 1000 # number of iterations
    burn = 0 # burn in 
    n_sim = n_iter - burn # number of stored simulations
    Z = array(0,c(n_sim,n)) # storage for cluster assignments
    PI = array(0,c(n_sim,K0)) # storage for cluster proportions
    PM = array(0,c(n_sim,K0,K0)) # storage for community connection params
    draw_logf_z_given_A = rep(0,n_sim) # P(z|Data) for MAP estimate of z
    
    for (i in 1:n_iter){
        
        # Step 1. pi
        ns = table(zs)
        as = ns + a0
        pis = rdirichlet(n=1, as)
        
        # Step 2. z
        for (i_sample in 1:n){ # loop through all nodes
            pi_star = rep(0,K0)
            for (k in 1:K0){ # loop through all cluster
                pi_star[k] = pis[k]
                for (j_sample in c(1:n)[-i_sample]){ # loop through all nodes except i_sample
                    pi_star[k] = pi_star[k] * 
                        (Ps[k,zs[j_sample]])^A[i_sample,j_sample] * 
                        (1.0-Ps[k,zs[j_sample]])^(1-A[i_sample,j_sample])
                } # for (j_sample in c(1:n_sample)[-i_sample])
            } # for (k in 1:K)
            pi_star = pi_star / sum(pi_star)
            zs[i_sample] = sample(x=(1:K0), size=1, prob=pi_star)
        } # for (i_sample in 1:n_sample)
        
        # Step 3. P
        for (a in 1:K0){
            for (b in a:K0){
                sum_A_ij = sum_one_minus_A_ij = 0 
                for (i_sample in 1:(n-1)){
                    for (j_sample in (i_sample+1):n){
                        if ( (zs[i_sample]==a)&(zs[j_sample]==b) ){ 
                            sum_A_ij = sum_A_ij + A[i_sample,j_sample]
                            sum_one_minus_A_ij = sum_one_minus_A_ij + (1-A[i_sample,j_sample])
                        } # if ( (z_vec[i_sample]==a)&(z_vec[j_sample]==b) )
                    } # for (j_sample in (i_sample+1):n_sample)
                } # for (i_sample in 1:n_sample)
                beta1_star = b10 + sum_A_ij ; beta2_star = b20 + sum_one_minus_A_ij
                Ps[a,b] = Ps[b,a] = rbeta(n=1,beta1_star,beta2_star) 
            } # for (b in (a+1):K)
        } # for (a in 1:(K-1))
        
        # Calculating P(z|Data) to find MAP
        logf_z_given_A = 0
        for (i_sample in 1:n){
            logf_z_given_A = logf_z_given_A + log(pis[zs[i_sample]])
        } # for (i_sample in 1:n_sample)
        for (i_sample in 1:(n-1)){
            for (j_sample in (i_sample+1):n){
                logf_z_given_A = logf_z_given_A + 
                    A[i_sample,j_sample]*log(Ps[zs[i_sample],zs[j_sample]]) + 
                    (1.0-A[i_sample,j_sample])*log(1.0-Ps[zs[i_sample],zs[j_sample]])
            } # for (j_sample in (i_sample+1):n_sample)
        } # for (j_sample in (i_sample+1):n_sample)
        
        # Store values
        if(i > burn)
        {
            j = i - burn
            Z[j,] = zs
            PI[j,] = pis
            PM[j,,] = Ps
            draw_logf_z_given_A[j] = logf_z_given_A
        }
        
        # Print status
        print(paste("Iteration",i))
        print(ns)
        if ( round(i/100)==(i/100) ){
            print( paste0("iter=",i) )
            print( "pi_vec=" )
            print( round(pis,3) )
        }
        
    }
    ret_list <- list(Z = Z,
                     PI = PI,
                     PM = PM,
                     logf = draw_logf_z_given_A,
                     z = Z[which.max(draw_logf_z_given_A),])
}