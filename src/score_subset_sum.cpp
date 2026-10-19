// src/score_subset_sum.cpp
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "compensated_arithmetic.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Sum g[subs[r, ]] for each row r (subs is 1-based from R).
// [[Rcpp::export]]
NumericVector rustat_subset_row_sums_g(IntegerMatrix subs, NumericVector g) {
  const int K = subs.nrow();
  const int s = subs.ncol();
  const int n = g.size();
  NumericVector out(K);
  int has_bad_index = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (K >= 256)
#endif
  for (int r = 0; r < K; ++r) {
    rdistnn_arithmetic::accumulator acc;
    bool valid = true;
    for (int c = 0; c < s; ++c) {
      const int index = subs(r, c);
      if (index == NA_INTEGER || index < 1 || index > n) {
#ifdef _OPENMP
#pragma omp critical(rustat_subset_row_sums_g_error)
#endif
        { has_bad_index = 1; }
        valid = false;
        break;
      }
      acc.add(g[index - 1]);
    }
    out[r] = valid ? acc.value() : NA_REAL;
  }
  if (has_bad_index) {
    Rcpp::stop("Index out of bounds in rustat_subset_row_sums_g.");
  }
  return out;
}
