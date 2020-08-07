#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sample_SBM_fast(NumericVector z, NumericMatrix P) 
{
    int n = z.size();
    NumericMatrix A(n);
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < i; j++)
        {
            int aij = R::rbinom(1,P(z[i] - 1,z[j] - 1));
            A(i,j) = aij;
            A(j,i) = aij;
        }
    }
    return A;
}
