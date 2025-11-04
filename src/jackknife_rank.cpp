// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// Helpers ---------------------------------------------------------------

inline int idx3(int t, int k, int q, int T, int S, int Q)
{
  // Flatten (t,k,q) with 0-based indices into [0, T*S*Q)
  // Layout: (((t)*S + k)*Q + q)
  return (((t)*S) + k) * Q + q;
}

// Fast delete-one jackknife for rank-weighted functionals (matrix phi)
// [[Rcpp::export]]
Rcpp::List rustat_jackknife_fast_rank(const IntegerMatrix &ranks,
                                      const NumericVector &a,
                                      const NumericMatrix &PHI,
                                      const bool normalize = true)
{
  const int T = ranks.nrow();
  const int K = ranks.ncol(); // must be >= s+1
  const int s = a.size();
  if (K < s + 1)
    stop("ncol(ranks) must be at least s+1.");
  const int n_ref = PHI.nrow();
  const int q = PHI.ncol();

  // Z = sum(a) or 1
  double Z = 0.0;
  if (normalize)
  {
    for (int k = 0; k < s; ++k)
      Z += a[k];
    if (Z == 0.0)
      stop("sum(a) == 0 with normalize=TRUE.");
  }
  else
  {
    Z = 1.0;
  }

  // Precompute dk[k] = a_{k+1} - a_k for k=0..s-2
  std::vector<double> dk(std::max(1, s - 1), 0.0);
  for (int k = 0; k <= s - 2; ++k)
    dk[k] = a[k] - a[k + 1];

  // Storage:
  // y_tk: (T x s x q) neighbors 1..s
  // y_s1: (T x q)    neighbor s+1
  std::vector<double> y_tk((size_t)T * s * q, 0.0);
  NumericMatrix y_s1(T, q);

  // S_t_divZ: (T x q) baseline per-target sums divided by Z
  NumericMatrix S_t_divZ(T, q);

  // Build y_tk, y_s1, and S_t_divZ
  for (int t = 0; t < T; ++t)
  {
    // neighbors 1..s
    for (int k = 0; k < s; ++k)
    {
      int j = ranks(t, k) - 1; // 0-based
      if (j < 0 || j >= n_ref)
        stop("ranks indices out of range.");
      for (int qq = 0; qq < q; ++qq)
      {
        double v = PHI(j, qq);
        y_tk[idx3(t, k, qq, T, s, q)] = v;
        S_t_divZ(t, qq) += a[k] * v;
      }
    }
    // neighbor s+1
    int jtail = ranks(t, s) - 1;
    if (jtail < 0 || jtail >= n_ref)
      stop("ranks indices out of range.");
    for (int qq = 0; qq < q; ++qq)
    {
      y_s1(t, qq) = PHI(jtail, qq);
    }
    // divide baseline by Z
    for (int qq = 0; qq < q; ++qq)
      S_t_divZ(t, qq) /= Z;
  }

  // theta0 = mean_t S_t_divZ[t,]
  NumericVector theta0(q);
  for (int t = 0; t < T; ++t)
  {
    for (int qq = 0; qq < q; ++qq)
      theta0[qq] += S_t_divZ(t, qq);
  }
  for (int qq = 0; qq < q; ++qq)
    theta0[qq] /= (double)T;

  // Precompute pre_t[t][r] = sum_{k=r}^{s-1} (a_{k+1}-a_k) * y_{t,k+1}
  // We store only for r=0..s-2; r=s-1 has zero.
  std::vector<double> pre_t((size_t)T * std::max(1, s - 1) * q, 0.0);
  if (s >= 2)
  {
    for (int t = 0; t < T; ++t)
    {
      // backward accumulation over k
      // acc = sum_{k=q}^{s-1} dk[k] * y_{t,k+1}
      // We'll fill for each r=q (0..s-2)
      for (int r = s - 2; r >= 0; --r)
      {
        // pre_t[t, r] = pre_t[t, r+1] + dk[r] * y_{t, r+1}
        for (int qq = 0; qq < q; ++qq)
        {
          double prev = (r == s - 2) ? 0.0 : pre_t[(size_t)t * (s - 1) * q + (r + 1) * q + qq];
          double add = dk[r] * y_tk[idx3(t, r + 1, qq, T, s, q)];
          pre_t[(size_t)t * (s - 1) * q + r * q + qq] = prev + add;
        }
      }
    }
  }

  // Build adjacency i -> list of (t,r) where r in 1..s+1
  std::vector<std::vector<std::pair<int, int>>> adj(n_ref);
  for (int t = 0; t < T; ++t)
  {
    for (int k = 0; k < s + 1; ++k)
    {
      int i = ranks(t, k) - 1;
      if (i < 0 || i >= n_ref)
        stop("ranks indices out of range.");
      adj[i].emplace_back(t, k + 1); // store r = 1..s+1
    }
  }

  // Output: theta_loo (n_ref x q), initialize to theta0
  NumericMatrix theta_loo(n_ref, q);
  for (int i = 0; i < n_ref; ++i)
  {
    for (int qq = 0; qq < q; ++qq)
      theta_loo(i, qq) = theta0[qq];
  }

  // Main loop over i: accumulate delta_sum and write theta_loo[i,]
  // Each i is independent -> parallelizable if OpenMP is available.
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n_ref; ++i)
  {
    std::vector<double> delta(q, 0.0);
    const auto &lst = adj[i];
    for (size_t m = 0; m < lst.size(); ++m)
    {
      int t = lst[m].first;
      int r = lst[m].second; // 1..s+1
      if (r > s)
        continue; // only affects if within support
      // base = S_t_divZ[t,]
      for (int qq = 0; qq < q; ++qq)
      {
        double base = S_t_divZ(t, qq);
        double term1 = (a[r - 1] / Z) * y_tk[idx3(t, r - 1, qq, T, s, q)];
        double term2 = 0.0;
        if (s >= 2 && r <= s - 1)
        {
          term2 = pre_t[(size_t)t * (s - 1) * q + (r - 1) * q + qq] / Z;
        }
        double term3 = (a[s - 1] / Z) * y_s1(t, qq);
        double Snew = base - term1 + term2 + term3;
        delta[qq] += (Snew - base);
      }
    }
    for (int qq = 0; qq < q; ++qq)
    {
      theta_loo(i, qq) = theta0[qq] + delta[qq] / (double)T;
    }
  }

  return Rcpp::List::create(
      _["theta0"] = theta0,
      _["theta_loo"] = theta_loo,
      _["info"] = Rcpp::List::create(
          _["s"] = s, _["T"] = T, _["q"] = q, _["normalize"] = normalize));
}
