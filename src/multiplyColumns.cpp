
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat multiplyColumns (arma::mat &A, arma::rowvec x) {
  arma::mat out = A;

  out.each_row() %= x;

  return out;
}
