// src/score_plqr.cpp
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Scores-only context for PLQR (no geometry, no localizer, no weights)
struct PlqrCtx
{
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

// [[Rcpp::export]]
SEXP rustat_plqr_ctx_create(NumericVector Y,
                            NumericMatrix D,
                            NumericVector qx,
                            NumericMatrix mD,
                            NumericVector beta,
                            double tau,
                            int j)
{
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

  XPtr<PlqrCtx> ctx(new PlqrCtx(), true);
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

// (Optional) update beta in-place without rebuilding the context
// [[Rcpp::export]]
void rustat_plqr_ctx_set_beta(SEXP ctx_, NumericVector beta)
{
  XPtr<PlqrCtx> ctx(ctx_);
  if (beta.size() != ctx->p_d)
    stop("beta length must equal ncol(D).");
  ctx->beta = beta;
}

// Compute all per-row PLQR contributions g (parallel)
// g_i = (tau - 1{ Y_i <= D_i^T beta + qx_i }) * (D_{i,j} - mD_{i,j})
// [[Rcpp::export]]
NumericVector rustat_plqr_g_all(SEXP ctx_)
{
  XPtr<PlqrCtx> ctx(ctx_);
  const int n = ctx->n;
  const int p_d = ctx->p_d;
  const int j0 = ctx->j - 1;
  NumericVector g(n);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n >= 1024)
#endif
  for (int i = 0; i < n; ++i)
  {
    double db = 0.0;
    for (int pd = 0; pd < p_d; ++pd)
      db += ctx->D(i, pd) * ctx->beta[pd];
    const double thr = db + ctx->qx[i];
    const int ind = (ctx->Y[i] <= thr) ? 1 : 0;
    const double resid = ctx->tau - ind;
    const double rDj = ctx->D(i, j0) - ctx->mD(i, j0);
    g[i] = resid * rDj;
  }
  return g;
}

// Same, but for a supplied 1-based index set (useful for block work)
// [[Rcpp::export]]
NumericVector rustat_plqr_g_indexed(IntegerVector idxs, SEXP ctx_)
{
  XPtr<PlqrCtx> ctx(ctx_);
  const int m = idxs.size();
  const int p_d = ctx->p_d;
  const int j0 = ctx->j - 1;
  NumericVector g(m);
  int first_bad_index = -1;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (m >= 1024)
#endif
  for (int t = 0; t < m; ++t)
  {
    const int i = idxs[t] - 1;
    if (i < 0 || i >= ctx->n) {
#ifdef _OPENMP
#pragma omp critical(rustat_plqr_g_indexed_error)
#endif
      {
        if (first_bad_index < 0) first_bad_index = i + 1;
      }
      g[t] = NA_REAL;
      continue;
    }
    double db = 0.0;
    for (int pd = 0; pd < p_d; ++pd)
      db += ctx->D(i, pd) * ctx->beta[pd];
    const double thr = db + ctx->qx[i];
    const int ind = (ctx->Y[i] <= thr) ? 1 : 0;
    const double resid = ctx->tau - ind;
    const double rDj = ctx->D(i, j0) - ctx->mD(i, j0);
    g[t] = resid * rDj;
  }
  if (first_bad_index >= 0) {
    Rcpp::stop("Index out of bounds in plqr_g_indexed.");
  }
  return g;
}
