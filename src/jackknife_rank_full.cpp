// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
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
// This is the algorithm described in the implementation note
// (markdown_notes/archive/completed_plans_2026-07-09/tdnn_rdistnn_implementation_plan.md, "Phase 2") and applies equally to
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

  // ---- theta0 : q x m -------------------------------------------------------
  NumericMatrix theta0(q, m);
  for (int t = 0; t < m; ++t) {
    for (int i = 0; i < n; ++i) {
      int r0 = ranks(i, t) - 1; // 0-based rank
      double wt = wn[r0];
      for (int qq = 0; qq < q; ++qq)
        theta0(qq, t) += wt * PHI(i, qq);
    }
  }

  // ---- theta_loo : dim c(n, q, m), column-major [i + n*qq + n*q*t] ---------
  NumericVector loo((R_xlen_t)n * q * m);

  // Per-target working buffers (reused across targets and signal dims)
  std::vector<int> ord(n);
  std::vector<double> P(n); // P[r0]: prefix sum for 0-based rank r0
  std::vector<double> S(n); // S[r0]: shifted suffix sum for 0-based rank r0

  for (int t = 0; t < m; ++t) {

    // Build ord: ord[r0] = 0-based obs index at 0-based rank r0.
    // Sort observation indices by their rank in column t.
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return ranks(a, t) < ranks(b, t); });

    for (int qq = 0; qq < q; ++qq) {

      // ---- forward prefix pass ----------------------------------------------
      // P[r0] = sum_{j=0}^{r0-1} wnm1[j] * z[j]   (0-based j, wnm1 size n-1)
      //
      // P[0] = 0
      // P[r0] = P[r0-1] + wnm1[r0-1] * z[r0-1]   for r0 = 1 .. n-1
      //
      // (P[n-1] is the last needed value; wnm1[r0-1] is valid for r0 <= n-1.)
      P[0] = 0.0;
      for (int r0 = 1; r0 < n; ++r0)
        P[r0] = P[r0 - 1] + wnm1[r0 - 1] * PHI(ord[r0 - 1], qq);

      // ---- backward shifted-suffix pass -------------------------------------
      // S[r0] = sum_{j=r0+1}^{n-1} wnm1[j-1] * z[j]
      //
      // S[n-1] = 0
      // S[r0]  = wnm1[r0] * z[r0+1] + S[r0+1]   for r0 = n-2 .. 0
      //
      // (wnm1[r0] is valid for r0 <= n-2.)
      S[n - 1] = 0.0;
      for (int r0 = n - 2; r0 >= 0; --r0)
        S[r0] = wnm1[r0] * PHI(ord[r0 + 1], qq) + S[r0 + 1];

      // ---- fill theta_loo ---------------------------------------------------
      // For obs i with 0-based rank ri = ranks(i,t) - 1:
      //   theta_loo[i, qq, t] = P[ri] + S[ri]
      R_xlen_t base = (R_xlen_t)qq * n + (R_xlen_t)qq * 0 /* see below */;
      // Correct flat index: i + n*qq + n*q*t
      R_xlen_t offset = (R_xlen_t)n * qq + (R_xlen_t)n * q * t;
      for (int i = 0; i < n; ++i) {
        int ri = ranks(i, t) - 1; // 0-based rank
        loo[offset + i] = P[ri] + S[ri];
      }
      (void)base; // suppress unused-variable warning
    }
  }

  loo.attr("dim") = IntegerVector::create(n, q, m);

  return List::create(_["theta0"] = theta0, _["theta_loo"] = loo,
                      _["info"] = List::create(_["n"] = n, _["m"] = m, _["q"] = q));
}
