// src/dist_l1_scaled.cpp
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <numeric>
#include <algorithm>
using namespace Rcpp;

// Core: row-wise scaled L1 distances
static inline void dist_l1_scaled_core(const NumericMatrix &X,
                                       const NumericVector &x0,
                                       const NumericVector &inv_s,
                                       std::vector<double> &out)
{
  const int n = X.nrow(), p = X.ncol();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n >= 1024)
#endif
  for (int i = 0; i < n; ++i)
  {
    double acc = 0.0;
    for (int j = 0; j < p; ++j)
    {
      const double d = (X(i, j) - x0[j]) * inv_s[j];
      acc += d >= 0.0 ? d : -d; // |d|
    }
    out[i] = acc;
  }
}

// [[Rcpp::export]]
NumericVector rustat_dist_l1_scaled(NumericMatrix X,
                                    NumericVector x0,
                                    NumericVector inv_s)
{
  if (x0.size() != X.ncol())
    stop("x0 length must equal ncol(X).");
  if (inv_s.size() != X.ncol())
    stop("inv_s length must equal ncol(X).");
  const int n = X.nrow();
  std::vector<double> buf(n);
  dist_l1_scaled_core(X, x0, inv_s, buf);
  return NumericVector(buf.begin(), buf.end());
}

// [[Rcpp::export]]
NumericMatrix rustat_dist_l1_scaled_batch(NumericMatrix X,
                                          NumericMatrix X0,
                                          NumericVector inv_s)
{
  if (X0.ncol() != X.ncol())
    stop("ncol(X0) must match ncol(X).");
  if (inv_s.size() != X.ncol())
    stop("inv_s length must equal ncol(X).");
  const int n = X.nrow(), m = X0.nrow();
  NumericMatrix out(n, m);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (m >= 8)
#endif
  for (int j = 0; j < m; ++j)
  {
    std::vector<double> col(n);
    dist_l1_scaled_core(X, X0.row(j), inv_s, col);
    for (int i = 0; i < n; ++i)
      out(i, j) = col[i];
  }
  return out;
}
