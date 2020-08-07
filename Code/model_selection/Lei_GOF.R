# Implement Lei's goodness of fit test given a cluster assignment vector
library(RMTstat)
Lei_GOF <- function(A,z,alpha = 0.05)
{
    n = length(z)
    Ns = table(z) # Number in each cluster 
    K0 = length(unique(z)) # putative number of clusters
    B <- matrix(0,nrow = K0,ncol = K0) # empty B matrix
    ids <- 1:n # node id vector
    for(k in 1:K0)
    {
        for(l in k:K0)
        {
            # inter-community connections
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
            # intra-community connections
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
    
    T_crit = qtw(1-(alpha/2))
    
    ret <- c(Tz,T_crit)
    names(ret) <- c("Test Statistic","Critical Value")
    return(ret)
}
