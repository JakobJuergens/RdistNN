// src/dist_l1.cpp
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "distance_arithmetic.hpp"

#include "parallel_threads.h"

#include <cmath>   // std::fabs
#include <cstddef> // std::size_t
#include <limits>

using namespace Rcpp;

// =====================================================================
// Core: row-wise L1 with optional per-dimension scaling.
// X is n x p (column-major), x0 length p. invsp may be nullptr.
// Writes into out[0..n-1].
// =====================================================================
static inline void dist_l1_core_ptr(const double *Xp, int n, int p, const double *x0p,
                                    const double *invsp, // may be nullptr
                                    double *out) {
  const std::size_t stride_n = static_cast<std::size_t>(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(rustat::thread_limit()) schedule(static) if (n >= 1024)
#endif
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    bool exceptional = false;
    for (int j = 0; j < p; ++j) {
      const double diff = Xp[i + j * stride_n] - x0p[j];
      const double d =
          invsp ? rustat::scaled_distance_difference(Xp[i + j * stride_n], x0p[j], invsp[j]) : diff;
      exceptional = exceptional || (invsp && (std::fabs(diff) >= rustat::distance_range_limit ||
                                              (diff != 0.0 && invsp[j] != 0.0 &&
                                               std::fabs(d) < std::numeric_limits<double>::min())));
      acc += std::fabs(d);
    }
    out[i] = exceptional ? rustat::scaled_norm_row(Xp, n, p, i, x0p, 1, invsp, false) : acc;
  }
}

// =====================================================================
// Batch helper: X is n x p, X0 is m x p. Writes n x m (column-major).
// invsp may be nullptr. Writes directly into 'out' (REAL of R matrix).
// =====================================================================
static inline void dist_l1_batch_impl(const double *Xp, int n, int p, const double *X0p, int m,
                                      const double *invsp, double *out) {
  const std::size_t stride_n = static_cast<std::size_t>(n);
  const std::size_t stride_m = static_cast<std::size_t>(m);

#ifdef _OPENMP
#pragma omp parallel for num_threads(rustat::thread_limit()) schedule(static) if (m >= 8)
#endif
  for (int q = 0; q < m; ++q) {
    for (int i = 0; i < n; ++i) {
      double acc = 0.0;
      bool exceptional = false;
      for (int j = 0; j < p; ++j) {
        const double x0j = X0p[q + j * stride_m];       // column-major X0
        const double diff = Xp[i + j * stride_n] - x0j; // column-major X
        const double d =
            invsp ? rustat::scaled_distance_difference(Xp[i + j * stride_n], x0j, invsp[j]) : diff;
        exceptional =
            exceptional || (invsp && (std::fabs(diff) >= rustat::distance_range_limit ||
                                      (diff != 0.0 && invsp[j] != 0.0 &&
                                       std::fabs(d) < std::numeric_limits<double>::min())));
        acc += std::fabs(d);
      }
      out[i + q * stride_n] =
          exceptional ? rustat::scaled_norm_row(Xp, n, p, i, X0p + q, m, invsp, false) : acc;
    }
  }
}

// =====================================================================
// Exports: unscaled/scaled single & batch
// =====================================================================

// [[Rcpp::export(rng = false)]]
NumericVector rustat_dist_l1(NumericMatrix X, NumericVector x0) {
  const int n = X.nrow(), p = X.ncol();
  if (x0.size() != p)
    stop("x0 length must equal ncol(X).");

  NumericVector out(n);
  dist_l1_core_ptr(REAL(X), n, p, REAL(x0), /*invsp=*/nullptr, REAL(out));
  return out;
}

// [[Rcpp::export(rng = false)]]
NumericVector rustat_dist_l1_scaled(NumericMatrix X, NumericVector x0, NumericVector inv_s) {
  const int n = X.nrow(), p = X.ncol();
  if (x0.size() != p)
    stop("x0 length must equal ncol(X).");
  if (inv_s.size() != p)
    stop("inv_s length must equal ncol(X).");

  NumericVector out(n);
  dist_l1_core_ptr(REAL(X), n, p, REAL(x0), REAL(inv_s), REAL(out));
  return out;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix rustat_dist_l1_batch(NumericMatrix X, NumericMatrix X0) {
  const int n = X.nrow(), p = X.ncol(), m = X0.nrow();
  if (X0.ncol() != p)
    stop("ncol(X0) must match ncol(X).");

  NumericMatrix out(n, m);
  dist_l1_batch_impl(REAL(X), n, p, REAL(X0), m, /*invsp=*/nullptr, REAL(out));
  return out;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix rustat_dist_l1_scaled_batch(NumericMatrix X, NumericMatrix X0, NumericVector inv_s) {
  const int n = X.nrow(), p = X.ncol(), m = X0.nrow();
  if (X0.ncol() != p)
    stop("ncol(X0) must match ncol(X).");
  if (inv_s.size() != p)
    stop("inv_s length must equal ncol(X).");

  NumericMatrix out(n, m);
  dist_l1_batch_impl(REAL(X), n, p, REAL(X0), m, REAL(inv_s), REAL(out));
  return out;
}
