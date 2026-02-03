// src/fast_order.cpp
// [[Rcpp::plugins(cpp17)]]
#include "fast_order.hpp"
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>

namespace rustat
{

  bool all_finite_nonneg(const double *x, int n) noexcept
  {
    for (int i = 0; i < n; ++i)
    {
      double v = x[i];
      if (!(std::isfinite(v)) || v < 0.0)
        return false;
    }
    return true;
  }

  struct StableCmp
  {
    const double *x;
    explicit StableCmp(const double *px) : x(px) {}

    // Rank among non-finites: -Inf(0) < +Inf(1) < NaN(2)
    static inline int nf_rank(double v)
    {
      if (std::isnan(v))
        return 2;
      // here v is non-finite and not NaN => +/-Inf
      return (v > 0.0) ? 1 : 0; // +Inf -> 1, -Inf -> 0
    }

    inline bool operator()(int i, int j) const
    {
      const double a = x[i], b = x[j];
      const bool ai = std::isfinite(a), aj = std::isfinite(b);

      // 1) Both finite -> numeric ascending, then index for stability
      if (ai && aj)
      {
        if (a < b)
          return true;
        if (b < a)
          return false;
        return i < j;
      }

      // 2) One finite, one non-finite -> finite first
      if (ai != aj)
        return ai; // TRUE if i is finite and j is not

      // 3) Both non-finite -> -Inf < +Inf < NaN, then index
      const int ra = nf_rank(a), rb = nf_rank(b);
      if (ra != rb)
        return ra < rb;
      return i < j; // deterministic tie-break
    }
  };

  static void argsort_radix64_nonneg(const double *x, int n, int *out_idx)
  {
    if (n <= 1)
    {
      if (n == 1)
        out_idx[0] = 0;
      return;
    }

    // Bulk bit copy: doubles -> uint64
    std::vector<std::uint64_t> keys(n);
    std::memcpy(keys.data(), x, sizeof(double) * (size_t)n);

    // Init permutation [0..n-1] and scratch buffer
    for (int i = 0; i < n; ++i)
      out_idx[i] = i;
    std::vector<int> buf(n);
    int *from = out_idx;
    int *to = buf.data();

    // 11-bit buckets => 6 passes cover 64 bits
    static constexpr int PASS_BITS = 11;
    static constexpr int RAD = 1 << PASS_BITS;                      // 2048
    static constexpr int PASSES = (64 + PASS_BITS - 1) / PASS_BITS; // 6
    const std::uint64_t MASK = (std::uint64_t)(RAD - 1);

    // heap-allocated counts to avoid large stack frames
    std::vector<size_t> counts(RAD);

    for (int pass = 0; pass < PASSES; ++pass)
    {
      std::fill(counts.begin(), counts.end(), 0);
      const int shift = pass * PASS_BITS;

      // count
      for (int i = 0; i < n; ++i)
      {
        const std::uint64_t k = keys[from[i]];
        const unsigned int b = (unsigned int)((k >> shift) & MASK);
        ++counts[b];
      }

      // prefix sums
      size_t sum = 0;
      for (int b = 0; b < RAD; ++b)
      {
        const size_t c = counts[b];
        counts[b] = sum;
        sum += c;
      }

      // stable scatter
      for (int i = 0; i < n; ++i)
      {
        const std::uint64_t k = keys[from[i]];
        const unsigned int b = (unsigned int)((k >> shift) & MASK);
        to[counts[b]++] = from[i];
      }

      std::swap(from, to);
    }

    if (from != out_idx)
    {
      std::memcpy(out_idx, from, sizeof(int) * (size_t)n);
    }
  }

  void argsort_indices(const double *x, int n, int *out_idx)
  {
    if (all_finite_nonneg(x, n))
    {
      argsort_radix64_nonneg(x, n, out_idx);
    }
    else
    {
      for (int i = 0; i < n; ++i)
        out_idx[i] = i;
      std::stable_sort(out_idx, out_idx + n, StableCmp(x));
    }
  }

  void ranks_from_values(const double *x, int n, int *ranks_out)
  {
    std::vector<int> ord(n);
    argsort_indices(x, n, ord.data()); // 0-based permutation
    for (int pos = 0; pos < n; ++pos)
    {
      int orig = ord[pos];
      ranks_out[orig] = pos + 1; // 1-based ranks
    }
  }

  int topk_indices(const double *x, int n, int k, int *out_idx_sorted)
  {
    if (k < 0)
      k = 0;
    if (k > n)
      k = n;
    if (k == 0)
      return 0;

    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i)
      idx[i] = i;
    auto cmp = StableCmp(x);

    std::nth_element(idx.begin(), idx.begin() + k, idx.end(), cmp);
    std::sort(idx.begin(), idx.begin() + k, cmp);

    for (int i = 0; i < k; ++i)
      out_idx_sorted[i] = idx[i];
    return k;
  }
} // namespace rustat

// ------------------ R wrappers (optional exports) ---------------------------

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector rustat_argsort_dist(NumericVector d)
{
  const int n = d.size();
  IntegerVector ord(n);
  if (n == 0)
    return ord;
  std::vector<int> idx(n);
  rustat::argsort_indices(REAL(d), n, idx.data());
  for (int i = 0; i < n; ++i)
    ord[i] = idx[i] + 1;
  return ord;
}

// [[Rcpp::export]]
IntegerVector rustat_ranks_from_dist(NumericVector d)
{
  const int n = d.size();
  IntegerVector ranks(n);
  if (n == 0)
    return ranks;
  std::vector<int> r(n);
  rustat::ranks_from_values(REAL(d), n, r.data());
  for (int i = 0; i < n; ++i)
    ranks[i] = r[i];
  return ranks;
}

// [[Rcpp::export]]
IntegerVector rustat_topk_dist(NumericVector d, int k)
{
  const int n = d.size();
  if (k < 0)
    k = 0;
  if (k > n)
    k = n;
  IntegerVector out(k);
  if (k == 0)
    return out;
  std::vector<int> idx(k);
  rustat::topk_indices(REAL(d), n, k, idx.data());
  for (int i = 0; i < k; ++i)
    out[i] = idx[i] + 1;
  return out;
}
