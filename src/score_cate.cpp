// src/score_cate.cpp
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "parallel_threads.h"
using namespace Rcpp;

// Context for CATE: per-row contribution g_i = rho_i
struct CateCtx {
  NumericVector rho; // length n (pseudo-outcome)
  int n;
};

// Create a context external pointer (weights/ordering handled elsewhere)
// [[Rcpp::export]]
SEXP rustat_cate_ctx_create(NumericVector rho) {
  if (rho.size() <= 0)
    stop("rho must be non-empty.");
  XPtr<CateCtx> ctx(new CateCtx(), true);
  ctx->rho = rho;
  ctx->n = rho.size();
  return ctx;
}

// Return ALL per-row contributions g (parallel-friendly read only)
// For CATE, g == rho, so we can just return the underlying vector.
// [[Rcpp::export]]
NumericVector rustat_cate_g_all(SEXP ctx_) {
  XPtr<CateCtx> ctx(ctx_);
  // Returning the SEXP is fine; it's read-only on C++ side.
  return ctx->rho;
}

// Return g on a set of indices (1-based from R) — filled in parallel
// Useful if you later want to do block processing.
// [[Rcpp::export]]
NumericVector rustat_cate_g_indexed(IntegerVector idxs, SEXP ctx_) {
  XPtr<CateCtx> ctx(ctx_);
  const int m = idxs.size();
  NumericVector out(m);
  int has_bad_index = 0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(rustat::thread_limit()) schedule(static) if (m >= 1024)
#endif
  for (int t = 0; t < m; ++t) {
    const int index = idxs[t];
    if (index == NA_INTEGER || index < 1 || index > ctx->n) {
#ifdef _OPENMP
#pragma omp critical(rustat_cate_g_indexed_error)
#endif
      { has_bad_index = 1; }
      out[t] = NA_REAL;
      continue;
    }
    const int i = index - 1; // 1-based -> 0-based
    out[t] = ctx->rho[i];
  }
  if (has_bad_index) {
    Rcpp::stop("Index out of bounds in cate_g_indexed.");
  }
  return out;
}
