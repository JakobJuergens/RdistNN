#ifndef RDISTNN_PARALLEL_THREADS_H
#define RDISTNN_PARALLEL_THREADS_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace rustat {

// Evaluate on the calling R thread, before entering each parallel region.
// num_threads() limits our own teams without altering another package's runtime.
inline int thread_limit() {
  SEXP option = Rf_GetOption1(Rf_install("RdistNN.num_threads"));
  int requested = 2;
  if (option != R_NilValue) {
    if ((TYPEOF(option) != INTSXP && TYPEOF(option) != REALSXP) || Rf_xlength(option) != 1)
      Rcpp::stop("`RdistNN.num_threads` must be a positive integer scalar.");
    const double value = Rcpp::as<double>(option);
    if (!std::isfinite(value) || value < 1 || value > std::numeric_limits<int>::max() ||
        value != std::floor(value))
      Rcpp::stop("`RdistNN.num_threads` must be a positive integer scalar.");
    requested = static_cast<int>(value);
  }
#ifdef _OPENMP
  if (omp_in_parallel())
    return 1;
  return std::min(requested, std::min(omp_get_max_threads(), omp_get_thread_limit()));
#else
  return 1;
#endif
}

} // namespace rustat
#endif
