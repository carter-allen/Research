# Model selection script for SBMs.

library(igraph)
library(sbmlogit)
library(sbmlhelpers)
library(RMTstat)

# simulate data from a K = 4 SBM with following connectivity matrix
n <- 200
p <- 0.75
q <- 0.05
K = 4
P <- matrix(c(p,q,q,q,
              q,p,q,q,
              q,q,p,q,
              q,q,q,p),
            nrow = 4, 
            ncol = 4,
            byrow = TRUE)
G1 <- sample_sbm(n = n,
                 block.sizes = rep(n/K,K),
                 pref.matrix = P,
                 directed = FALSE,
                 loops = FALSE)
G1_g <- as_adjacency_matrix(G1,
                          type = "both",
                          sparse = FALSE)

# implement SVT to recover K
S <- svd(G1_g)
s <- S$d
sum(s > sqrt(n))

# implement Lei's goodness of fit test to recover K
K0 <- 4
fitK <- sbmlogit.mcmc(G1,
                      alpha = K0,
                      nsamples = 1000)
z <- get_labels(fitK)
Ns <- table(z)
A <- G1_g
B <- matrix(0,nrow = K0,ncol = K0)
ids <- 1:n
for(k in 1:K0)
{
    for(l in k:K0)
    {
        if(k != l)
        {
            nk = Ns[k]
            nl = Ns[l]
            is = ids[z == k]
            js = ids[z == l]
            b = 0
            for(i in is)
            {
                for(j in js)
                {
                    b = b + A[i,j]
                }
            }
            B[k,l] = B[l,k] = b/(nk*nl)
        }
        else
        {
            nk = Ns[k]
            is = ids[z == k]
            b = 0
            for(i in 1:length(is))
            {
                for(j in i:length(is)) 
                {
                    b = b + A[is[i],is[j]]
                }
            }
            B[k,l] = B[l,k] = b/((nk*(nk-1))/2)
        }
    }
}

Az = A
for(i in 1:ncol(A))
{
    for(j in 1:nrow(A))
    {
        if(i == j)
        {
            A[i,j] = 0
        }
        else
        {
            Az[i,j] = (A[i,j] - B[z[i],z[j]])/(sqrt((n-1)*B[z[i],z[j]]*(1-B[z[i],z[j]])))
        }
    }
}

Sz = svd(Az)
sz = max(Sz$d)
Tz = (n^(2/3))*(sz-2)

a = 0.05
T_crit = qtw(1-(a/2))
