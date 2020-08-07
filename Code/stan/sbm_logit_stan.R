# R code for trying to fit SBM logit in stan

# simulate data
n <- 100 # number of nodes
N <- choose(n,2) # number of possible edges
K <- 3 # number of communities
pi <- rep(1/K,K) # marginal community probs
z <- sample(x = c(1,2,3),
            size = n,
            replace = TRUE,
            prob = pi) # generate cluster node counts 
Z <- make_Z_mat(z) # construct community association matrix
Z_lower <- Z[lower.tri(Z)] # construct long form associations

# construct community probability matrix
P <- matrix(c(0.50,0.05,0.05,
              0.05,0.50,0.05,
              0.05,0.05,0.50),
            nrow = 3,
            ncol = 3,
            byrow = TRUE)

# get probability of edge for each possible pair of nodes
theta <- rep(0,N)
for(i in 1:N)
{
    zij = Z_lower[i]
    zi = as.numeric(substr(zij,1,1))
    zj = as.numeric(substr(zij,2,2))
    theta[i] = P[zi,zj]
}

# sample edges according to theta
y = rbinom(N,1,theta)

# construct design matrix
X_cols = sort(unique(Z_lower))
X <- matrix(0,nrow = N, ncol = length(unique(Z_lower)))
for(i in 1:nrow(X))
{
    for(j in 1:ncol(X))
    {
        X[i,j] = ifelse(Z_lower[i] == X_cols[j],
                        1,
                        0)
    }
}

s
