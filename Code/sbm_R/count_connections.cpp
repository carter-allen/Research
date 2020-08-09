#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void count_connections(NumericVector zs, NumericMatrix A) {
  int K0 = unique(zs).length();
  int n = zs.length();
  NumericMatrix C(K0);
  for(int i = 0; i < n; i++)
  {
      for(int j = 0; j < i; j++)
      {
          C(zs[i]-1,zs[j]-1) += A(i,j);
      }
  }
  Rcout << C << std::endl;
}

