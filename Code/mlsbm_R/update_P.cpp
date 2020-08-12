#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix update_P(List A,
                       NumericVector zs,
                       int K0,
                       double b10,
                       double b20) 
{
  int L = A.length();
  int n = zs.size();
  NumericMatrix P(K0);
  for(int a = 0; a < K0; a ++)
  {
    for(int b = a; b < K0; b ++)
    {
      int sum_A_ij = 0;
      int sum_one_minus_A_ij = 0;
      for(int i_sample = 0; i_sample < n-1; i_sample ++)
      {
        for(int j_sample = i_sample + 1; j_sample < n; j_sample ++)
        {
          if(zs[i_sample] == (a+1) && zs[j_sample] == (b+1))
          {
            for(int l = 0; l < L; l ++)
            {
              NumericMatrix Al = A[l];
              sum_A_ij = sum_A_ij + Al(i_sample,j_sample);
              sum_one_minus_A_ij = sum_one_minus_A_ij + (1-Al(i_sample,j_sample));
            }
          }
        }
      }
      double beta1_star = b10 + sum_A_ij ; 
      double beta2_star = b20 + sum_one_minus_A_ij;
      //Rcout << a << "," << b << " : " <<sum_A_ij << "; " << sum_one_minus_A_ij << std::endl;
      NumericVector Pab = Rcpp::rbeta(1,beta1_star,beta2_star);
      double pab = Pab(0);
      P(a,b) = pab;
      P(b,a) = pab;
    }
  }
  return P;
}
