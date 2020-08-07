// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void update_z(NumericVector z_s, 
                       NumericMatrix Amat, 
                       NumericVector pi_s, 
                       NumericMatrix P_s,
                       int K0,
                       NumericVector classes) 
{
  int n = z_s.size();
  NumericVector z_ret = z_s;
  for(int i_sample = 0; i_sample < n; i_sample ++)
  {
      NumericVector pi_star (K0);
      for(int k = 0; k < K0; k ++)
      {
          pi_star[k] = pi_s[k];
          for(int j_sample = 0; j_sample < n; j_sample ++)
          {
              if(i_sample != j_sample)
              {
                  //Rcout  << pi_star[k] << ": " << pow(P_s(k,z_s[j_sample]),Amat(i_sample,j_sample)) << ": " << pow((1 - P_s(k,z_s[j_sample])),(1-Amat(i_sample,j_sample))) <<  std::endl;
                  pi_star[k] = pi_star[k] * 
                      pow(P_s(k,z_s[j_sample]),Amat(i_sample,j_sample)) * 
                      pow((1 - P_s(k,z_s[j_sample])),(1-Amat(i_sample,j_sample)));
              }
          }
      }
      Rcout << i_sample << ":" << pi_star << std::endl;
      pi_star = pi_star / sum(pi_star);
      //Rcout << i_sample << ": " << z_s[i_sample] << ": " << pi_star << std::endl;
      NumericVector z = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pi_star);
      z_ret[i_sample] = z[0];
      //Rcout << i_sample << ": " << z << std::endl;
  }
}

