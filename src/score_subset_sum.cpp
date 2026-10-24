// src/score_subset_sum.cpp
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "compensated_arithmetic.hpp"
#include "parallel_threads.h"
using namespace Rcpp;

// Sum g[subs[r, ]] for each row r (subs is 1-based from R).
// [[Rcpp::export(rng = false)]]
NumericVector rustat_subset_row_sums_g(IntegerMatrix subs, NumericVector g) {
  const int K = subs.nrow();
  const int s = subs.ncol();
  const int n = g.size();
  NumericVector out(K);
  int has_bad_index = 0;
  // Materialize R buffers on the calling thread before entering workers.
  const int *subsets = INTEGER(subs);
  const double *contributions = REAL(g);
  double *sums = REAL(out);
  const double missing = NA_REAL;

#ifdef _OPENMP
#pragma omp parallel for num_threads(rustat::thread_limit()) schedule(static) if (K >= 256)
#endif
  for (int r = 0; r < K; ++r) {
    rdistnn_arithmetic::accumulator acc;
    bool valid = true;
    for (int c = 0; c < s; ++c) {
      const int index = subsets[r + static_cast<R_xlen_t>(K) * c];
      if (index == NA_INTEGER || index < 1 || index > n) {
#ifdef _OPENMP
#pragma omp critical(rustat_subset_row_sums_g_error)
#endif
        { has_bad_index = 1; }
        valid = false;
        break;
      }
      acc.add(contributions[index - 1]);
    }
    sums[r] = valid ? acc.value() : missing;
  }
  if (has_bad_index) {
    Rcpp::stop("Index out of bounds in rustat_subset_row_sums_g.");
  }
  return out;
}
