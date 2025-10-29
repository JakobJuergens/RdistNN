// src/dist_l1.cpp
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

// Core: row-wise L1 distances to x0
static inline void dist_l1_core(const NumericMatrix &X,
                                const NumericVector &x0,
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
      const double dif = X(i, j) - x0[j];
      acc += dif >= 0.0 ? dif : -dif; // fabs without function call
    }
    out[i] = acc;
  }
}

// [[Rcpp::export]]
NumericVector rustat_dist_l1(NumericMatrix X, NumericVector x0)
{
  if (x0.size() != X.ncol())
    stop("x0 length must equal ncol(X).");
  const int n = X.nrow();
  std::vector<double> buf(n);
  dist_l1_core(X, x0, buf);
  return NumericVector(buf.begin(), buf.end());
}

// [[Rcpp::export]]
List rustat_order_by_l1(NumericMatrix X, NumericVector x0)
{
  if (x0.size() != X.ncol())
    stop("x0 length must equal ncol(X).");
  const int n = X.nrow();
  std::vector<double> dx(n);
  dist_l1_core(X, x0, dx);

  // order by (dx, index) — deterministic ties like base R
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b)
            {
              if (dx[a] < dx[b])
                return true;
              if (dx[a] > dx[b])
                return false;
              return a < b; // tie-break by index
            });

  IntegerVector ord1(n);
  for (int i = 0; i < n; ++i)
    ord1[i] = ord[i] + 1;
  return List::create(_["dx"] = NumericVector(dx.begin(), dx.end()),
                      _["ord"] = ord1);
}
