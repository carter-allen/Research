# function to contruct the community association matrix for 
# a given community labeling vector z.
# If z is nx1 then Z is nxn where each element (i,j) represents 
# the communities of nodes i and j.

make_Z_mat <- function(z)
{
    n = length(z)
    Z = matrix(0,nrow = n, ncol = n)
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            zi = z[i]
            zj = z[j]
            zs <- sort(c(zi,zj))
            Z[i,j] = paste0(zs[1],zs[2])
        }
    }
    return(Z)
}