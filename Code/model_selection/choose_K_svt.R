# Implement SVT to choose K given symmetric adjacency matrix
choose_K_svt <- function(A)
{
    if(!(isSymmetric(A) & is.numeric(A)))
    {
        return("Please give a symmetric adjacency matrix")
    }
    else
    {
        n <- nrow(A)
        S <- svd(A)
        s <- S$d
        K <- sum(s > sqrt(n))
        return(K)
    }
}