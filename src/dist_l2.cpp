// src/dist_l2.cpp
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "fast_order.hpp"
#include "distance_arithmetic.hpp"

#include "parallel_threads.h"

#include <vector>
#include <cmath>   // std::sqrt
#include <cstddef> // std::size_t
#include <limits>

using namespace Rcpp;

// Common intermediate boundary: Stata reserves doubles >= 2^1023 as missing.
static constexpr double l2_range_limit = rustat::distance_range_limit;

// Round each squared coordinate before adding it, as in the R fallback.
// The volatile boundary prevents fused multiply-add from changing exact ties.
static inline double rounded_square(const double d) {
  volatile double product = d * d;
  return product;
}

// Exceptional recomputation only: scalar workspace, no per-row buffer.
static inline double l2_range_row(const double *Xp, int n, int p, int i, const double *x0p,
                                  int stride0, const double *invsp, bool squared) {
  if (invsp)
    return rustat::scaled_norm_row(Xp, n, p, i, x0p, stride0, invsp, true, squared);
  double scale = 0.0, ssq = 1.0;
  for (int j = 0; j < p; ++j) {
    const double x = Xp[i + static_cast<std::size_t>(j) * n];
    const double target = x0p[static_cast<std::size_t>(j) * stride0];
    double d = x - target;
    const double a = std::fabs(d);
    if (!std::isfinite(a))
      return a;
    if (a != 0.0) {
      if (scale < a) {
        const double ratio = scale / a;
        const volatile double term = ssq * rounded_square(ratio);
        ssq = 1.0 + term;
        scale = a;
      } else
        ssq += rounded_square(a / scale);
    }
  }
  return squared ? (scale * ssq) * scale : scale * std::sqrt(ssq);
}

// =====================================================================
// Single core for row-wise L2 distances with optional per-dim scaling.
// If invsp == nullptr we do unscaled L2; otherwise scaled by inv_s[j].
// X is n x p (column-major), x0 is length p. Results into out[0..n-1].
// =====================================================================
static inline void dist_l2_core_ptr(const double *Xp, int n, int p, const double *x0p,
                                    const double *invsp, // may be nullptr
                                    double *out, const bool squared) {
  const std::size_t stride_n = static_cast<std::size_t>(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(rustat::thread_limit()) schedule(static) if (n >= 1024)
#endif
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    bool exceptional = false;
    for (int j = 0; j < p; ++j) {
      const double diff = Xp[i + j * stride_n] - x0p[j];
      const double d = invsp ? diff * invsp[j] : diff;
      exceptional = exceptional || (invsp && std::fabs(diff) >= l2_range_limit);
      exceptional = exceptional || (invsp && diff != 0.0 && invsp[j] != 0.0 &&
                                    std::fabs(d) < std::numeric_limits<double>::min());
      const double square = rounded_square(d);
      // Nonfinite squares propagate to acc; check that once after the loop.
      exceptional = exceptional || (d != 0.0 && square < std::numeric_limits<double>::min());
      acc += square;
    }
    out[i] = (exceptional || !std::isfinite(acc) || acc >= l2_range_limit)
                 ? l2_range_row(Xp, n, p, i, x0p, 1, invsp, squared)
                 : (squared ? acc : std::sqrt(acc));
  }
}

// =====================================================================
// Batch helper: X is n x p, X0 is m x p. Writes an n x m result matrix.
// invsp may be nullptr for unscaled. (We write directly into 'out'.)
// =====================================================================
static inline void dist_l2_batch_impl(const double *Xp, int n, int p, const double *X0p, int m,
                                      const double *invsp, double *out, bool squared) {
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
        const double d = invsp ? diff * invsp[j] : diff;
        exceptional = exceptional || (invsp && std::fabs(diff) >= l2_range_limit);
        exceptional = exceptional || (invsp && diff != 0.0 && invsp[j] != 0.0 &&
                                      std::fabs(d) < std::numeric_limits<double>::min());
        const double square = rounded_square(d);
        // Nonfinite squares propagate to acc; check that once after the loop.
        exceptional = exceptional || (d != 0.0 && square < std::numeric_limits<double>::min());
        acc += square;
      }
      out[i + q * stride_n] = (exceptional || !std::isfinite(acc) || acc >= l2_range_limit)
                                  ? l2_range_row(Xp, n, p, i, X0p + q, m, invsp, squared)
                                  : (squared ? acc : std::sqrt(acc));
    }
  }
}

// =====================================================================
// Thin wrappers: un/scaled, single/batch, and order-by distance
// =====================================================================

// [[Rcpp::export(rng = false)]]
NumericVector rustat_dist_l2(NumericMatrix X, NumericVector x0, bool squared = false) {
  const int n = X.nrow(), p = X.ncol();
  if (x0.size() != p)
    stop("x0 length must equal ncol(X).");

  NumericVector out(n);
  dist_l2_core_ptr(REAL(X), n, p, REAL(x0), /*invsp=*/nullptr, out.begin(), squared);
  return out;
}

// [[Rcpp::export(rng = false)]]
NumericVector rustat_dist_l2_scaled(NumericMatrix X, NumericVector x0, NumericVector inv_s,
                                    bool squared = false) {
  const int n = X.nrow(), p = X.ncol();
  if (x0.size() != p)
    stop("x0 length must equal ncol(X).");
  if (inv_s.size() != p)
    stop("inv_s length must equal ncol(X).");

  NumericVector out(n);
  dist_l2_core_ptr(REAL(X), n, p, REAL(x0), REAL(inv_s), out.begin(), squared);
  return out;
}

// Optional convenience: unscaled batch
// [[Rcpp::export(rng = false)]]
NumericMatrix rustat_dist_l2_batch(NumericMatrix X, NumericMatrix X0, bool squared = false) {
  const int n = X.nrow(), p = X.ncol(), m = X0.nrow();
  if (X0.ncol() != p)
    stop("ncol(X0) must match ncol(X).");

  NumericMatrix out(n, m);
  dist_l2_batch_impl(REAL(X), n, p, REAL(X0), m,
                     /*invsp=*/nullptr, REAL(out), squared);
  return out;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix rustat_dist_l2_scaled_batch(NumericMatrix X, NumericMatrix X0, NumericVector inv_s,
                                          bool squared = false) {
  const int n = X.nrow(), p = X.ncol(), m = X0.nrow();
  if (X0.ncol() != p)
    stop("ncol(X0) must match ncol(X).");
  if (inv_s.size() != p)
    stop("inv_s length must equal ncol(X).");

  NumericMatrix out(n, m);
  dist_l2_batch_impl(REAL(X), n, p, REAL(X0), m, REAL(inv_s), REAL(out), squared);
  return out;
}

// Shared ordering helper
static inline IntegerVector order_by_dx(const std::vector<double> &dx) {
  const int n = static_cast<int>(dx.size());
  std::vector<int> ord0(n);
  rustat::argsort_indices(dx.data(), n, ord0.data());
  IntegerVector ord(n);
  for (int i = 0; i < n; ++i)
    ord[i] = ord0[i] + 1; // 1-based for R
  return ord;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List rustat_order_by_l2(NumericMatrix X, NumericVector x0, bool squared = false) {
  const int n = X.nrow(), p = X.ncol();
  if (x0.size() != p)
    stop("x0 length must equal ncol(X).");

  std::vector<double> dx(n);
  dist_l2_core_ptr(REAL(X), n, p, REAL(x0), /*invsp=*/nullptr, dx.data(), squared);
  return List::create(_["dx"] = NumericVector(dx.begin(), dx.end()), _["ord"] = order_by_dx(dx));
}

// [[Rcpp::export(rng = false)]]
Rcpp::List rustat_order_by_l2_scaled(NumericMatrix X, NumericVector x0, NumericVector inv_s,
                                     bool squared = false) {
  const int n = X.nrow(), p = X.ncol();
  if (x0.size() != p)
    stop("x0 length must equal ncol(X).");
  if (inv_s.size() != p)
    stop("inv_s length must equal ncol(X).");

  std::vector<double> dx(n);
  dist_l2_core_ptr(REAL(X), n, p, REAL(x0), REAL(inv_s), dx.data(), squared);
  return List::create(_["dx"] = NumericVector(dx.begin(), dx.end()), _["ord"] = order_by_dx(dx));
}
