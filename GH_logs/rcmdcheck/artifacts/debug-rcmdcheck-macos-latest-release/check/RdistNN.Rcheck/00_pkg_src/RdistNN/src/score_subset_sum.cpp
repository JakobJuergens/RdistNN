// src/score_subset_sum.cpp
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
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
  int first_bad_index = -1;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (K >= 256)
#endif
  for (int r = 0; r < K; ++r) {
    double acc = 0.0;
    for (int c = 0; c < s; ++c) {
      const int i = subs(r, c) - 1;
      if (i < 0 || i >= n) {
#ifdef _OPENMP
#pragma omp critical(rustat_subset_row_sums_g_error)
#endif
        {
          if (first_bad_index < 0)
            first_bad_index = i + 1;
        }
        acc = NA_REAL;
        break;
      }
      acc += g[i];
    }
    out[r] = acc;
  }
  if (first_bad_index >= 0) {
    Rcpp::stop("Index out of bounds in rustat_subset_row_sums_g.");
  }
  return out;
}
