count_connections_R <- function(zs,A)
{
    n = length(zs)
    K0 = length(unique(zs))
    C = matrix(0,K0,K0)
    for(i in 1:n)
    {
        for(j in 1:i)
        {
            C[z[i],z[j]] = C[z[i],z[j]] + A[i,j]
        }
    }
    return(C)
}