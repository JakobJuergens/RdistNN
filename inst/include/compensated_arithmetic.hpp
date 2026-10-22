#ifndef RDISTNN_COMPENSATED_ARITHMETIC_HPP
#define RDISTNN_COMPENSATED_ARITHMETIC_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// Mirrored in the Stata companion. Finite double-precision arithmetic, with
// separately rounded products; no wider-type or parallel-reduction policy.
namespace rdistnn_arithmetic {

struct accumulator {
  double high = 0.0;
  double low = 0.0;

  void add(double value) {
    const double next = high + value;
    if (!std::isfinite(next)) {
      high = next;
      low = 0.0;
      return;
    }
    low += std::abs(high) >= std::abs(value) ? (high - next) + value : (value - next) + high;
    high = next;
  }

  void add_product(double x, double y) {
    volatile double product = x * y;
    add(product);
  }

  double value() const {
    return high + low;
  }
};

inline double sum(const double *x, int n) {
  accumulator total;
  for (int i = 0; i < n; ++i)
    total.add(x[i]);
  return total.value();
}

inline double dot(const double *x, const double *y, int n) {
  accumulator total;
  for (int i = 0; i < n; ++i)
    total.add_product(x[i], y[i]);
  return total.value();
}

// Sorted nonnegative weights. Return a 1-based index, or zero for invalid mass.
// Keeping both scans here makes boundary selection independent of long double.
inline int weighted_cdf_index(const double *weights, int n, double prob) {
  const double mass = sum(weights, n);
  if (!std::isfinite(mass) || mass <= 0)
    return 0;
  accumulator prefix;
  for (int i = 0; i < n; ++i) {
    prefix.add(weights[i]);
    if (prefix.value() / mass >= prob)
      return i + 1;
  }
  return n;
}

// z is rank-ordered and support-centered; table has n-1 entries. Reuse two
// n-element arrays: forward high/low state, then overwrite high with results.
// The reverse state stays scalar, so retaining corrections adds no n-array.
inline void deletion_sums(const double *z, const double *table, int n, double offset, double *high,
                          double *low) {
  if (n == 0)
    return;
  high[0] = low[0] = 0.0;
  accumulator left;
  for (int r = 1; r < n; ++r) {
    left.add_product(table[r - 1], z[r - 1]);
    high[r] = left.high;
    low[r] = left.low;
  }
  accumulator right;
  for (int r = n; r-- > 0;) {
    accumulator merged;
    merged.add(high[r]);
    merged.add(right.high);
    merged.add(low[r]);
    merged.add(right.low);
    merged.add(offset);
    high[r] = merged.value();
    if (r > 0)
      right.add_product(table[r - 1], z[r]);
  }
}

// Canonical n x (m*q) output, in target/component blocks. Callers validate
// permutations, shapes and finite inputs before entering this computation.
inline void rank_loo(const int *ranks, const double *signal, const double *table, int n, int m,
                     int q, double *out) {
  const double mass = sum(table, n - 1);
  double weight_amplitude = 0.0;
  for (int a = 0; a < n - 1; ++a)
    weight_amplitude = std::max(weight_amplitude, std::abs(table[a]));
  std::vector<int> order(n);
  std::vector<double> high(n), low(n), z(n);
  for (int t = 0; t < m; ++t) {
    const int *r = ranks + static_cast<std::size_t>(t) * n;
    for (int i = 0; i < n; ++i)
      order[r[i] - 1] = i;
    for (int component = 0; component < q; ++component) {
      const double *values = signal + static_cast<std::size_t>(component) * n;
      double minimum = std::numeric_limits<double>::infinity();
      double maximum = -std::numeric_limits<double>::infinity();
      int first = n;
      for (int a = 0; a < n; ++a) {
        const bool active = (a < n - 1 && table[a] != 0) || (a > 0 && table[a - 1] != 0);
        if (active) {
          if (first == n)
            first = order[a];
          minimum = std::min(minimum, values[order[a]]);
          maximum = std::max(maximum, values[order[a]]);
        }
      }
      double level = first == n ? 0.0 : values[first];
      const double limit = std::ldexp(1.0, 1023);
      if (first < n && !(std::abs(maximum - level) < limit && std::abs(minimum - level) < limit))
        level = minimum > 0 ? minimum : (maximum < 0 ? maximum : 0.0);
      // Match the shared R/Mata underflow-amplitude and nonzero-mass policy.
      const double amplitude = std::max(std::abs(minimum), std::abs(maximum));
      if (first < n &&
          (amplitude < std::numeric_limits<double>::min() ||
           weight_amplitude * amplitude < std::numeric_limits<double>::min() ||
           (mass != 0.0 && (std::abs(mass) < std::numeric_limits<double>::min() ||
                            std::abs(mass) * amplitude < std::numeric_limits<double>::min()))))
        level = 0.0;
      for (int a = 0; a < n; ++a) {
        const bool active = (a < n - 1 && table[a] != 0) || (a > 0 && table[a - 1] != 0);
        if (!active)
          continue;
        const double value = values[order[a]];
        const double difference = value - level;
        const double back = difference - value;
        const double error = (value - (difference - back)) - (level + back);
        if (error != 0) {
          level = 0;
          break;
        }
      }
      for (int a = 0; a < n; ++a) {
        const bool active = (a < n - 1 && table[a] != 0) || (a > 0 && table[a - 1] != 0);
        z[a] = active ? values[order[a]] - level : 0.0;
      }
      volatile double offset = level * mass;
      deletion_sums(z.data(), table, n, offset, high.data(), low.data());
      double *destination = out + (static_cast<std::size_t>(t) * q + component) * n;
      for (int a = 0; a < n; ++a)
        destination[order[a]] = high[a];
    }
  }
}

// Average actual rank weights, never signals. The padded endpoints are only
// used for blocks that cannot occur before/after the deletion, respectively.
inline double block_table_mean(const double *table, int n, int first, int last) {
  const int count = last - first;
  if (count == 0)
    return 0.0;
  accumulator total;
  for (int r = first; r < last; ++r)
    total.add(r >= 0 && r < n - 1 ? table[r] : 0.0);
  const double mass = total.value();
  if (std::isfinite(mass))
    return mass / count;
  accumulator divided;
  for (int r = first; r < last; ++r)
    divided.add(r >= 0 && r < n - 1 ? table[r] / count : 0.0);
  return divided.value();
}

// Direct actual-observation deletions for noncanonical supplied groups, or
// the requested reference algorithm. This uncommon supplied-rank boundary is
// O(n^2*q); ordinary covariate-aware ranks use the linear block algorithm.
inline void block_loo_direct(const std::vector<int> &order, const int *groups, const double *signal,
                             const double *table, int n, int q, double *out) {
  std::vector<int> counts(n), rank_of_observation(n);
  std::vector<accumulator> masses(n);
  std::vector<double> weights(n);
  for (int r = 0; r < n; ++r)
    rank_of_observation[order[r]] = r;
  for (int deleted = 0; deleted < n; ++deleted) {
    std::fill(counts.begin(), counts.end(), 0);
    std::fill(masses.begin(), masses.end(), accumulator{});
    // Match public duplicate-weight averaging: sum each block's assigned
    // table entries in original observation order, even for supplied ranks.
    for (int observation = 0; observation < n; ++observation) {
      const int r = rank_of_observation[observation];
      if (r == deleted)
        continue;
      const int group = groups[observation] - 1;
      ++counts[group];
      masses[group].add(table[r - (r > deleted)]);
    }
    bool needs_divided = false;
    for (int group = 0; group < n; ++group) {
      weights[group] = counts[group] ? masses[group].value() / counts[group] : 0.0;
      needs_divided = needs_divided || !std::isfinite(weights[group]);
    }
    if (needs_divided) {
      std::fill(masses.begin(), masses.end(), accumulator{});
      for (int observation = 0; observation < n; ++observation) {
        const int r = rank_of_observation[observation];
        const int group = groups[observation] - 1;
        if (r != deleted && !std::isfinite(weights[group]))
          masses[group].add(table[r - (r > deleted)] / counts[group]);
      }
      for (int group = 0; group < n; ++group)
        if (!std::isfinite(weights[group]))
          weights[group] = masses[group].value();
    }
    for (int component = 0; component < q; ++component) {
      accumulator total;
      for (int r = 0; r < n; ++r)
        if (r != deleted)
          total.add_product(signal[order[r] + static_cast<std::size_t>(component) * n],
                            weights[groups[order[r]] - 1]);
      out[order[deleted] + static_cast<std::size_t>(component) * n] = total.value();
    }
  }
}

// Centered ORIGINAL signals, not rounded block means. Ranks/groups are 1-based;
// output has column-major target/component blocks, as rank_loo(). Canonical
// contiguous groups take O(n*m*q) work and O(n) scratch. Prefix/suffix states
// strictly outside the deleted block combine with its actual survivor states;
// no large removed-row contribution is subtracted from a rounded block mean.
inline void block_rank_loo(const int *ranks, const double *signal, const double *table,
                           const int *groups, int n, int m, int q, double *out,
                           bool direct = false) {
  for (int i = 0; i < n; ++i)
    if (groups[i] < 1 || groups[i] > n)
      throw std::invalid_argument("groups must contain integer identifiers in 1..n");
  std::vector<int> order(n), starts, ends, seen(n);
  std::vector<double> before, after, survivor;
  for (int target = 0; target < m; ++target) {
    const int *ranks_t = ranks + static_cast<std::size_t>(target) * n;
    std::fill(seen.begin(), seen.end(), 0);
    for (int i = 0; i < n; ++i) {
      const int represented_rank = ranks_t[i];
      if (represented_rank < 1 || represented_rank > n)
        throw std::invalid_argument("ranks must be permutations of 1..n");
      const int rank = represented_rank - 1;
      if (seen[rank]++)
        throw std::invalid_argument("ranks must be permutations of 1..n");
      order[rank] = i;
    }
    starts.clear();
    ends.clear();
    std::fill(seen.begin(), seen.end(), 0);
    bool canonical = true;
    for (int a = 0; a < n;) {
      const int group = groups[order[a]] - 1;
      if (seen[group]++)
        canonical = false;
      int b = a + 1;
      while (b < n && groups[order[b]] - 1 == group) {
        if (order[b] < order[b - 1])
          canonical = false;
        ++b;
      }
      starts.push_back(a);
      ends.push_back(b);
      a = b;
    }
    double *target_out = out + static_cast<std::size_t>(target) * n * q;
    if (direct || !canonical) {
      block_loo_direct(order, groups, signal, table, n, q, target_out);
      continue;
    }
    const int blocks = static_cast<int>(starts.size());
    before.resize(blocks);
    after.resize(blocks);
    survivor.resize(blocks);
    for (int g = 0; g < blocks; ++g) {
      before[g] = block_table_mean(table, n, starts[g], ends[g]);
      after[g] = block_table_mean(table, n, starts[g] - 1, ends[g] - 1);
      survivor[g] = block_table_mean(table, n, starts[g], ends[g] - 1);
    }
    std::vector<accumulator> prefix(blocks + 1), suffix(blocks + 1), local(n);
    for (int component = 0; component < q; ++component) {
      const double *z = signal + static_cast<std::size_t>(component) * n;
      prefix[0] = suffix[blocks] = accumulator{};
      for (int g = 0; g < blocks; ++g) {
        prefix[g + 1] = prefix[g];
        for (int r = starts[g]; r < ends[g]; ++r)
          prefix[g + 1].add_product(z[order[r]], before[g]);
      }
      for (int g = blocks; g-- > 0;) {
        suffix[g] = suffix[g + 1];
        for (int r = ends[g]; r-- > starts[g];)
          suffix[g].add_product(z[order[r]], after[g]);
      }
      for (int g = 0; g < blocks; ++g) {
        accumulator left, right;
        for (int r = starts[g]; r < ends[g]; ++r) {
          local[r] = left;
          left.add_product(z[order[r]], survivor[g]);
        }
        for (int r = ends[g]; r-- > starts[g];) {
          accumulator total;
          total.add(prefix[g].high);
          total.add(suffix[g + 1].high);
          total.add(local[r].high);
          total.add(right.high);
          total.add(prefix[g].low);
          total.add(suffix[g + 1].low);
          total.add(local[r].low);
          total.add(right.low);
          target_out[order[r] + static_cast<std::size_t>(component) * n] = total.value();
          right.add_product(z[order[r]], survivor[g]);
        }
      }
    }
  }
}

} // namespace rdistnn_arithmetic
#endif
