// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "compensated_arithmetic.hpp"
using namespace Rcpp;

// Generic O(n * m * q) rank-order delete-one jackknife via prefix/suffix sums.
//
// For any rank-weight rule fully specified by:
//   wn[r]   (1-based r = 1..n):   weight for rank r at full sample size n,
//   wnm1[r] (1-based r = 1..n-1): weight for rank r at leave-one-out size n-1,
//
// the delete-one estimator for target t deleting obs i with 1-based rank r_i is:
//
//   theta_loo[i, t] = P[r_i] + S[r_i]
//
// where (using 1-based rank-order positions and 0-based C++ arrays):
//   z[r] = PHI(ord[r], qq)          signal at rank-r observation
//
//   P[r] = sum_{j=1}^{r-1} wnm1[j] * z[j]
//        = cumulative prefix at wnm1-scale up to (but not including) r
//
//   S[r] = sum_{j=r+1}^{n} wnm1[j-1] * z[j]
//        = shifted suffix (after deletion rank r, obs at rank j > r moves to j-1)
//
// Both passes are O(n) per (t, qq), giving O(n * m * q) total.
//
// The rank-table deletion identity is documented in
// markdown_notes/design/14_integrated_tuning_grid_reuse.md and applies to
// exact DNN, asymptotic DNN, exact TDNN, and asymptotic TDNN — any rule
// whose rank-weight table can be precomputed outside this function.
//
// Arguments
// ---------
// ranks : n x m integer matrix (1-based ranks, each column a permutation of 1..n)
// PHI   : n x q numeric signal matrix
// wn    : length-n numeric vector, rank-weight table at sample size n
// wnm1  : length-(n-1) numeric vector, rank-weight table at sample size n-1
//
// Return value
// ------------
// A named List with:
//   theta0    : q x m NumericMatrix, baseline estimates sum_i wn[r_i]*phi_i
//   theta_loo : NumericVector with dim attribute c(n, q, m),
//               element [i, qq, t] at position i + n*qq + n*q*t (0-based, C order)
//   info      : List(n, m, q)
//
// [[Rcpp::export]]
Rcpp::List rustat_jackknife_rank_full(const IntegerMatrix &ranks, const NumericMatrix &PHI,
                                      const NumericVector &wn, const NumericVector &wnm1) {
  const int n = ranks.nrow();
  const int m = ranks.ncol();
  const int q = PHI.ncol();

  if (PHI.nrow() != n)
    stop("nrow(PHI) must equal nrow(ranks).");
  if ((int)wn.size() != n)
    stop("length(wn) must equal nrow(ranks).");
  if ((int)wnm1.size() != n - 1)
    stop("length(wnm1) must equal nrow(ranks) - 1.");

  NumericMatrix theta0(q, m);
  for (int t = 0; t < m; ++t) {
    for (int component = 0; component < q; ++component) {
      rdistnn_arithmetic::accumulator total;
      for (int i = 0; i < n; ++i)
        total.add_product(wn[ranks(i, t) - 1], PHI(i, component));
      theta0(component, t) = total.value();
    }
  }
  NumericVector loo(static_cast<R_xlen_t>(n) * q * m);
  rdistnn_arithmetic::rank_loo(ranks.begin(), PHI.begin(), wnm1.begin(), n, m, q, loo.begin());
  loo.attr("dim") = IntegerVector::create(n, q, m);
  return List::create(_["theta0"] = theta0, _["theta_loo"] = loo,
                      _["info"] = List::create(_["n"] = n, _["m"] = m, _["q"] = q));
}
