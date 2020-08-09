print <- function(x,...)
{
    UseMethod("print",x)
}
print.SBM <- function(fit)
{
    n = length(fit$z)
    K = fit$K
    k = length(unique(fit$z))
    cat("Summary of Bayesian SBM fit","\n")
    cat("Number of specified clusters:",K,"\n")
    cat("Number of inferred clusters:",k,"\n")
    cat("Number of nodes in each cluster:","\n")
    cat(table(fit$z))
}