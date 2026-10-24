// src/score_plqr.cpp
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "parallel_threads.h"
#include "linear_threshold_comparison.hpp"
#include <cmath>
using namespace Rcpp;

// Scores-only context for PLQR (no geometry, no localizer, no weights)
struct PlqrCtx {
  NumericMatrix D;    // n x p_d
  NumericMatrix mD;   // n x p_d
  NumericVector Y;    // n
  NumericVector qx;   // n
  NumericVector beta; // p_d
  int n;              // rows
  int p_d;            // cols in D/mD
  int j;              // 1-based component index
  double tau;         // quantile in (0,1)
};

// [[Rcpp::export(rng = false)]]
SEXP rustat_plqr_ctx_create(NumericVector Y, NumericMatrix D, NumericVector qx, NumericMatrix mD,
                            NumericVector beta, double tau, int j) {
  const int n = Y.size();
  if (D.nrow() != n || mD.nrow() != n || qx.size() != n)
    stop("Y, D, mD, and qx must share the same number of rows n.");
  if (D.ncol() != mD.ncol())
    stop("D and mD must have the same number of columns.");
  if (beta.size() != D.ncol())
    stop("beta length must equal ncol(D).");
  if (!(j >= 1 && j <= D.ncol()))
    stop("j is out of range.");
  if (!(tau > 0.0 && tau < 1.0))
    stop("tau must be in (0,1).");

  XPtr<PlqrCtx> ctx(new PlqrCtx(), true, Rf_install("RdistNN::PlqrCtx"));
  ctx->D = D;
  ctx->mD = mD;
  ctx->Y = Y;
  ctx->qx = qx;
  ctx->beta = beta;
  ctx->n = n;
  ctx->p_d = D.ncol();
  ctx->j = j;
  ctx->tau = tau;
  return ctx;
}

// Compute all per-row PLQR contributions g (parallel)
// g_i = (tau - 1{ Y_i <= D_i^T beta + qx_i }) * (D_{i,j} - mD_{i,j})
// [[Rcpp::export(rng = false)]]
NumericVector rustat_plqr_g_all(SEXP ctx_) {
  // External pointers do not survive serialization. Check identity and address
  // before interpreting any payload, including a different real native pointer.
  if (TYPEOF(ctx_) != EXTPTRSXP || R_ExternalPtrTag(ctx_) != Rf_install("RdistNN::PlqrCtx") ||
      R_ExternalPtrAddr(ctx_) == nullptr)
    stop("Invalid PLQR context; recreate the score context before native evaluation.");
  const PlqrCtx *ctx = static_cast<PlqrCtx *>(R_ExternalPtrAddr(ctx_));
  const int n = ctx->n;
  const int p_d = ctx->p_d;
  const int j0 = ctx->j - 1;
  NumericVector g(n);
  // Resolve R objects on the calling thread; workers use only plain buffers.
  const double *design = REAL(ctx->D), *coefficients = REAL(ctx->beta);
  const double *nuisance = REAL(ctx->mD), *qx = REAL(ctx->qx), *outcome = REAL(ctx->Y);
  double *scores = REAL(g);
  const double tau = ctx->tau;
  const double missing = NA_REAL;

#ifdef _OPENMP
#pragma omp parallel for num_threads(rustat::thread_limit()) schedule(static) if (n >= 1024)
#endif
  for (int i = 0; i < n; ++i) {
    double db = 0.0;
    for (int pd = 0; pd < p_d; ++pd)
      db += design[i + static_cast<R_xlen_t>(n) * pd] * coefficients[pd];
    const double thr = db + qx[i];
    const int ind = rdistnn_arithmetic::linear_threshold_indicator(design + i, n, coefficients, p_d,
                                                                   qx[i], outcome[i], thr);
    if (ind < 0) {
      scores[i] = missing;
      continue;
    }
    const double resid = tau - ind;
    const double rDj =
        design[i + static_cast<R_xlen_t>(n) * j0] - nuisance[i + static_cast<R_xlen_t>(n) * j0];
    scores[i] = resid * rDj;
  }
  return g;
}
