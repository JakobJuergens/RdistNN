// src/dist_l2.cpp
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

// Core compute: L2 or squared L2 distances row-wise
static inline void dist_l2_core(const NumericMatrix& X,
                                const NumericVector& x0,
                                std::vector<double>& out,
                                const bool squared)
{
  const int n = X.nrow(), p = X.ncol();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n >= 1024)
#endif
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j) {
      const double zij = X(i, j) - x0[j];
      acc += zij * zij;
    }
    out[i] = squared ? acc : std::sqrt(acc);
  }
}

// [[Rcpp::export]]
NumericVector rustat_dist_l2(NumericMatrix X, NumericVector x0, bool squared = false) {
  if (x0.size() != X.ncol()) Rcpp::stop("x0 length must equal ncol(X).");
  const int n = X.nrow();
  std::vector<double> buf(n);
  dist_l2_core(X, x0, buf, squared);
  return NumericVector(buf.begin(), buf.end());
}

// [[Rcpp::export]]
Rcpp::List rustat_order_by_l2(NumericMatrix X, NumericVector x0, bool squared = false) {
  if (x0.size() != X.ncol()) Rcpp::stop("x0 length must equal ncol(X).");
  const int n = X.nrow();
  std::vector<double> dx(n);
  dist_l2_core(X, x0, dx, squared);

  // order by (distance, index) — deterministic ties (like R's order)
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){
    if (dx[a] < dx[b]) return true;
    if (dx[a] > dx[b]) return false;
    return a < b; // tie-break by index
  });

  IntegerVector ord1(n);
  for (int i = 0; i < n; ++i) ord1[i] = ord[i] + 1;
  return Rcpp::List::create(
    Rcpp::_["dx"]   = NumericVector(dx.begin(), dx.end()),
    Rcpp::_["ord"]  = ord1
  );
}
