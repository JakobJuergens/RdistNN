// src/dist_l2_scaled.cpp
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
#include <cmath>

using namespace Rcpp;

// ---------------------------------------------------------------------
// Core: row-wise scaled L2 or squared L2 (pointer-based, column-major)
// out[i] = || (X[i,] - x0) * inv_s ||_2  [or squared if 'squared']
// ---------------------------------------------------------------------
static inline void dist_l2_scaled_core_ptr(const double* Xp, int n, int p,
                                           const double* x0p,
                                           const double* invsp,
                                           std::vector<double>& out,
                                           const bool squared)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n >= 1024)
#endif
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j) {
      const double d = (Xp[i + j * n] - x0p[j]) * invsp[j];
      acc += d * d;
    }
    out[i] = squared ? acc : std::sqrt(acc);
  }
}

// ---------------------------------------------------------------------
// Exported: distances for a single x0
// ---------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector rustat_dist_l2_scaled(NumericMatrix X,
                                    NumericVector x0,
                                    NumericVector inv_s,
                                    bool squared = false)
{
  const int p = X.ncol();
  if (x0.size() != p)   stop("x0 length must equal ncol(X).");
  if (inv_s.size() != p) stop("inv_s length must equal ncol(X).");

  const int n = X.nrow();
  std::vector<double> dx(n);
  dist_l2_scaled_core_ptr(REAL(X), n, p, REAL(x0), REAL(inv_s), dx, squared);
  return NumericVector(dx.begin(), dx.end());
}

// ---------------------------------------------------------------------
// Exported: batch distances for X0 (m x p); no x0 copies per query
// ---------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix rustat_dist_l2_scaled_batch(NumericMatrix X,
                                          NumericMatrix X0,
                                          NumericVector inv_s,
                                          bool squared = false)
{
  const int n = X.nrow(), p = X.ncol(), m = X0.nrow();
  if (X0.ncol() != p)   stop("ncol(X0) must match ncol(X).");
  if (inv_s.size() != p) stop("inv_s length must equal ncol(X).");

  NumericMatrix out(n, m);
  const double* Xp   = REAL(X);
  const double* X0p  = REAL(X0);
  const double* invp = REAL(inv_s);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (m >= 8)
#endif
  for (int q = 0; q < m; ++q) {
    std::vector<double> col(n);
    // x0 component j is at X0p[q + j*m] (column-major)
    // Compute distances for this query q without copying the row.
    for (int i = 0; i < n; ++i) {
      double acc = 0.0;
      for (int j = 0; j < p; ++j) {
        const double x0j = X0p[q + j * (std::size_t)m];
        const double d   = (Xp[i + j * n] - x0j) * invp[j];
        acc += d * d;
      }
      out(i, q) = squared ? acc : std::sqrt(acc);
    }
  }
  return out;
}

// ---------------------------------------------------------------------
// Exported: order rows by scaled L2 (stable; ties by index)
// Returns both distances and the 1-based order
// ---------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rustat_order_by_l2_scaled(NumericMatrix X,
                                     NumericVector x0,
                                     NumericVector inv_s,
                                     bool squared = false)
{
  const int p = X.ncol();
  if (x0.size() != p)   stop("x0 length must equal ncol(X).");
  if (inv_s.size() != p) stop("inv_s length must equal ncol(X).");

  const int n = X.nrow();
  std::vector<double> dx(n);
  dist_l2_scaled_core_ptr(REAL(X), n, p, REAL(x0), REAL(inv_s), dx, squared);

  std::vector<int> ord0(n);
  rustat::argsort_indices(dx.data(), n, ord0.data());  // radix or stable fallback

  IntegerVector ord(n);
  for (int i = 0; i < n; ++i) ord[i] = ord0[i] + 1;

  return List::create(
    _["dx"]  = NumericVector(dx.begin(), dx.end()),
    _["ord"] = ord
  );
}
