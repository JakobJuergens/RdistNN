// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "compensated_arithmetic.hpp"
#include "parallel_threads.h"
using namespace Rcpp;

// Helpers ---------------------------------------------------------------

inline std::size_t idx3(int t, int k, int q, int T, int S, int Q) {
  // Flatten (t,k,q) with 0-based indices into [0, T*S*Q)
  // Layout: (((t)*S + k)*Q + q)
  return ((static_cast<std::size_t>(t) * S) + k) * Q + q;
}

// Fast delete-one jackknife for rank-weighted functionals (matrix phi)
// [[Rcpp::export]]
Rcpp::List rustat_jackknife_fast_rank(const IntegerMatrix &ranks, const NumericVector &a,
                                      const NumericMatrix &PHI, const bool normalize = true) {
  const int T = ranks.nrow();
  const int K = ranks.ncol(); // must be >= s+1
  const int s = a.size();
  if (K < s + 1)
    stop("ncol(ranks) must be at least s+1.");
  const int n_ref = PHI.nrow();
  const int q = PHI.ncol();

  const double mass = rdistnn_arithmetic::sum(a.begin(), s);
  const double Z = normalize ? mass : 1.0;
  if (normalize && Z == 0)
    stop("sum(a) == 0 with normalize=TRUE.");
  NumericMatrix baseline(T, q);
  std::vector<double> local(static_cast<std::size_t>(T) * (s + 1) * q);
  std::vector<rdistnn_arithmetic::accumulator> totals(q);
  std::vector<std::vector<std::pair<int, int>>> affected(n_ref);
  std::vector<int> local_ranks(s + 1);
  for (int k = 0; k <= s; ++k)
    local_ranks[k] = k + 1;
  std::vector<double> z(s + 1), values(s + 1);
  for (int t = 0; t < T; ++t) {
    for (int k = 0; k <= s; ++k) {
      const int index = ranks(t, k);
      if (index == NA_INTEGER || index < 1 || index > n_ref)
        stop("ranks indices out of range.");
      const int i = index - 1;
      if (k < s)
        affected[i].emplace_back(t, k);
    }
    for (int component = 0; component < q; ++component) {
      rdistnn_arithmetic::accumulator base;
      for (int k = 0; k <= s; ++k) {
        z[k] = PHI(ranks(t, k) - 1, component);
        if (k < s)
          base.add_product(a[k], z[k]);
      }
      baseline(t, component) = base.value() / Z;
      totals[component].add(baseline(t, component));
      rdistnn_arithmetic::rank_loo(local_ranks.data(), z.data(), a.begin(), s + 1, 1, 1,
                                   values.data());
      for (int k = 0; k <= s; ++k)
        local[idx3(t, k, component, T, s + 1, q)] = values[k] / Z;
    }
  }
  NumericVector theta0(q);
  NumericMatrix theta_loo(n_ref, q);
  for (int component = 0; component < q; ++component)
    theta0[component] = totals[component].value() / T;
#ifdef _OPENMP
#pragma omp parallel for num_threads(rustat::thread_limit())
#endif
  for (int i = 0; i < n_ref; ++i) {
    for (int component = 0; component < q; ++component) {
      rdistnn_arithmetic::accumulator total;
      total.add(totals[component].high);
      total.add(totals[component].low);
      for (const auto &location : affected[i]) {
        const int t = location.first, k = location.second;
        total.add(-baseline(t, component));
        total.add(local[idx3(t, k, component, T, s + 1, q)]);
      }
      theta_loo(i, component) = total.value() / T;
    }
  }
  return List::create(
      _["theta0"] = theta0, _["theta_loo"] = theta_loo,
      _["info"] = List::create(_["s"] = s, _["T"] = T, _["q"] = q, _["normalize"] = normalize));
}
