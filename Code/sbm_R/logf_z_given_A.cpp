#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double logf_z_given_A(NumericVector zs, NumericMatrix A, NumericVector pis, NumericMatrix Ps) {
  double logf_z_given_A = 0;
  int n = zs.length();
  for(int i = 0; i < n; i++)
  {
      logf_z_given_A = logf_z_given_A + log(pis[zs[i]-1]);
  }
  for(int i = 0; i < n-1; i++)
  {
      for(int j = i + 1; j < n; j++)
      {
          logf_z_given_A = logf_z_given_A + A(i,j)*log(Ps(zs[i]-1,zs[j]-1)) + (1.0-A(i,j))*log(1.0-Ps(zs[i]-1,zs[j]-1));
      }
  }
  return logf_z_given_A;
}