#include <Rcpp.h>
#include "rank_weight_tables.hpp"
using namespace Rcpp;
// [[Rcpp::export(rng = false)]]
NumericVector rustat_exact_rank_table(int n, int k, int s) {
  if (n < 1 || k < 1 || k > s || s > n)
    stop("Require 1 <= k <= s <= n.");
  const auto table = rdistnn_weights::exact_table(n, k, s);
  return wrap(table);
}
