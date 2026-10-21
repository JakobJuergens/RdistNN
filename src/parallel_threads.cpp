// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include "parallel_threads.h"

// Internal diagnostic: establish the actual native team and unchanged runtime
// limit in tests, including builds without OpenMP support.
// [[Rcpp::export]]
Rcpp::List rustat_thread_info() {
  const int limit = rustat::thread_limit();
  int team_size = 1;
  int runtime_max = 1;
  int runtime_limit = 1;
  bool openmp = false;
#ifdef _OPENMP
  openmp = true;
  runtime_max = omp_get_max_threads();
  runtime_limit = omp_get_thread_limit();
#pragma omp parallel num_threads(rustat::thread_limit())
  {
#pragma omp single
    team_size = omp_get_num_threads();
  }
#endif
  return Rcpp::List::create(Rcpp::_["openmp"] = openmp, Rcpp::_["limit"] = limit,
                            Rcpp::_["team_size"] = team_size, Rcpp::_["runtime_max"] = runtime_max,
                            Rcpp::_["runtime_limit"] = runtime_limit);
}
