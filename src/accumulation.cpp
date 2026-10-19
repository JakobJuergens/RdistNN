// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "compensated_arithmetic.hpp"

// [[Rcpp::export]]
double rustat_compensated_sum(const Rcpp::NumericVector &x) {
  return rdistnn_arithmetic::sum(x.begin(), x.size());
}

// [[Rcpp::export]]
int rustat_weighted_cdf_index(const Rcpp::NumericVector &w, double prob) {
  const int index = rdistnn_arithmetic::weighted_cdf_index(w.begin(), w.size(), prob);
  return index == 0 ? NA_INTEGER : index;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rustat_compensated_crossprod(const Rcpp::NumericMatrix &x,
                                                 const Rcpp::NumericMatrix &y) {
  if (x.nrow() != y.nrow())
    Rcpp::stop("Weighted signals require aligned rows.");
  const int n = x.nrow();
  Rcpp::NumericMatrix out(x.ncol(), y.ncol());
  for (int j = 0; j < y.ncol(); ++j)
    for (int k = 0; k < x.ncol(); ++k)
      out(k, j) = rdistnn_arithmetic::dot(x.begin() + static_cast<std::size_t>(k) * n,
                                          y.begin() + static_cast<std::size_t>(j) * n, n);
  return out;
}
