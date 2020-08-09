#' Implement SVT to choose K given symmetric adjacency matrix
#'
#' Function to estimate number of communities K from an adjacency matrix A
#' by computing the number of singular values above sqrt(n), where n is the number
#' of nodes in the graph.
#' @param A A symmetric adjacency matrix with entries of 0 or 1
#' @keywords SBM, SVT,
#' @export
#' @examples
#' choose_K_svt(A)
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