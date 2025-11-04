// src/dist_l1.cpp
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
#include <cmath> // std::fabs

using namespace Rcpp;

// ---------------------------------------------------------------------
// Core compute: L1 (Manhattan) distances row-wise (parallel-friendly)
// NOTE: Uses column-major indexing for speed + correctness.
// ---------------------------------------------------------------------
static inline void dist_l1_core(const NumericMatrix &X,
                                const NumericVector &x0,
                                std::vector<double> &out)
{
  const int n = X.nrow(), p = X.ncol();
  const double* Xp  = REAL(X);
  const double* x0p = REAL(x0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n >= 1024)
#endif
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    // column-major: element (i,j) is Xp[i + j*n]
    for (int j = 0; j < p; ++j) {
      acc += std::fabs(Xp[i + j * n] - x0p[j]);
    }
    out[i] = acc;
  }
}

// ---------------------------------------------------------------------
// Exported: just the distance vector (L1)
// ---------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector rustat_dist_l1(NumericMatrix X, NumericVector x0)
{
  if (x0.size() != X.ncol())
    stop("x0 length must equal ncol(X).");

  const int n = X.nrow();
  std::vector<double> dx(n);
  dist_l1_core(X, x0, dx);
  return NumericVector(dx.begin(), dx.end());
}

// ---------------------------------------------------------------------
// Exported: order rows by L1 distance (stable; ties by index)
// Returns both distances and the 1-based order
// ---------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rustat_order_by_l1(NumericMatrix X, NumericVector x0)
{
  if (x0.size() != X.ncol())
    stop("x0 length must equal ncol(X).");

  const int n = X.nrow();
  std::vector<double> dx(n);
  dist_l1_core(X, x0, dx);

  // Fast stable argsort (radix path for finite, nonnegative distances)
  std::vector<int> ord0(n);
  rustat::argsort_indices(dx.data(), n, ord0.data());  // 0-based

  IntegerVector ord(n);
  for (int i = 0; i < n; ++i) ord[i] = ord0[i] + 1;     // 1-based for R

  return List::create(
    _["dx"]  = NumericVector(dx.begin(), dx.end()),
    _["ord"] = ord
  );
}
