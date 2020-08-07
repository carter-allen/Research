// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp ;

// [[Rcpp::export]]
CharacterVector csample_char( CharacterVector x, 
                              int size,
                              bool replace, 
                              NumericVector prob = NumericVector::create()
) {
    CharacterVector ret = RcppArmadillo::sample(x, size, replace, prob) ;
    return ret ;
}