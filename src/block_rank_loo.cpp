// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "compensated_arithmetic.hpp"

// Internal original-signal block deletion kernel. R validates fitted state;
// this boundary also protects direct native callers before output writes.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector rustat_block_rank_loo(const Rcpp::IntegerMatrix &ranks,
                                          const Rcpp::NumericMatrix &signal,
                                          const Rcpp::NumericVector &table,
                                          const Rcpp::IntegerVector &groups, bool direct = false) {
  const int n = ranks.nrow(), m = ranks.ncol(), q = signal.ncol();
  if (n < 2 || m < 1 || q < 1 || signal.nrow() != n || table.size() != n - 1 || groups.size() != n)
    Rcpp::stop("Block deletions require aligned signals, groups, ranks and n-1 weights.");
  for (double value : signal)
    if (!std::isfinite(value))
      Rcpp::stop("Block deletion signals must be finite.");
  for (double value : table)
    if (!std::isfinite(value))
      Rcpp::stop("Block deletion weights must be finite.");
  Rcpp::NumericVector out(static_cast<R_xlen_t>(n) * m * q);
  rdistnn_arithmetic::block_rank_loo(ranks.begin(), signal.begin(), table.begin(), groups.begin(),
                                     n, m, q, out.begin(), direct);
  out.attr("dim") = Rcpp::IntegerVector::create(n, q, m);
  return out;
}
