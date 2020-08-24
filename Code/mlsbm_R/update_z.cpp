// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector update_z(NumericVector zs, 
                       List A, 
                       NumericMatrix Ps,
                       NumericVector pis,
                       NumericVector classes) 
{
    int L = A.length();
    int n = zs.length();
    int K0 = Ps.ncol();
    NumericVector z_ret(n);
    for(int i = 0; i < n; i ++)
    {
        NumericVector pi_star(K0);
        for(int k = 0; k < K0; k++)
        {
            pi_star[k] = log(pis[k]);
            for(int j = 0; j < n; j++)
            {
                if(i != j)
                {
                    for(int l = 0; l < L; l++)
                    {
                        NumericMatrix Al = A[l];
                        pi_star[k] = ((pi_star[k]) + 
                            log(pow(Ps(k,zs[j]-1),Al(i,j))) + 
                            log(pow((1 - Ps(k,zs[j]-1)),(1-Al(i,j)))));
                    }
                }
            }
            //Rcout << i << "," << k << ":" << pi_star[k] << std::endl;
        }
        pi_star = exp(pi_star - mean(pi_star)) / sum(exp(pi_star - mean(pi_star)));
        //Rcout << i << ":" << pi_star << std::endl;
        NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pi_star);
        z_ret[i] = zi[0];
    }
    return z_ret;
}


