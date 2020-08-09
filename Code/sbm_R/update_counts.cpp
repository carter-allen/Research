#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector update_counts(NumericVector zs, int K0) {
    NumericVector counts(K0);
    int n = zs.length();
    for(int k = 1; k <= K0; k++)
    {
        for(int i = 0; i < n; i++)
        {
            if(zs[i] == k)
            {
                counts[k-1] = counts[k-1] + 1;
            }
        }
    }
    return counts;
}
