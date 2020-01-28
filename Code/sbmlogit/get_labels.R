# function to obtain remapped posterior centroid estimator of cluster labels
# see Peng & Car. 2016
library(sbmlogit)

# most probable group label
mp = function(vec, K)
{
    v = rep(1:K)
    l = length(vec)
    
    for (i in 1:K)
    {
        v[i] = sum(vec==i)/l
    }
    return(v)
}

get_labels <- function(fit)
{
    K = fit$ngroups
    Sigma = fit$sample # posterior samples
    sigma = apply(t(apply(Sigma, 2, mp, K)), 1, which.max) # posterior estimator
    scentroid = sbmlogit.remap(sigma) # remapped posterior estimator
    return(scentroid)
}