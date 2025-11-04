// src/dist_l2.cpp
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "fast_order.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath> // std::sqrt

using namespace Rcpp;

// ---------------------------------------------------------------------
// Core compute: L2 or squared L2 distances row-wise (parallel-friendly)
// ---------------------------------------------------------------------
static inline void dist_l2_core(const NumericMatrix &X,
                                const NumericVector &x0,
                                std::vector<double> &out,
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

// ---------------------------------------------------------------------
// Exported: just the distance vector (L2 or squared L2)
// ---------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector rustat_dist_l2(NumericMatrix X, NumericVector x0, bool squared = false)
{
  if (x0.size() != X.ncol())
    stop("x0 length must equal ncol(X).");

  const int n = X.nrow();
  std::vector<double> dx(n);
  dist_l2_core(X, x0, dx, squared);
  return NumericVector(dx.begin(), dx.end());
}

// ---------------------------------------------------------------------
// Exported: order rows by L2 distance (stable; ties by index)
// Returns both distances and the 1-based order
// ---------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rustat_order_by_l2(NumericMatrix X, NumericVector x0, bool squared = false)
{
  if (x0.size() != X.ncol())
    stop("x0 length must equal ncol(X).");

  const int n = X.nrow();
  std::vector<double> dx(n);
  dist_l2_core(X, x0, dx, squared);

  // Fast stable argsort via rustat::argsort_indices (radix for nonneg finite)
  std::vector<int> ord0(n);
  rustat::argsort_indices(dx.data(), n, ord0.data());  // 0-based indices

  IntegerVector ord(n);
  for (int i = 0; i < n; ++i) ord[i] = ord0[i] + 1;     // 1-based for R

  return List::create(
    _["dx"]  = NumericVector(dx.begin(), dx.end()),
    _["ord"] = ord
  );
}
