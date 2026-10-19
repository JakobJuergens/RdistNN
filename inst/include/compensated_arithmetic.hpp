#ifndef RDISTNN_COMPENSATED_ARITHMETIC_HPP
#define RDISTNN_COMPENSATED_ARITHMETIC_HPP

#include <algorithm>
#include <cmath>
#include <limits>
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

} // namespace rdistnn_arithmetic
#endif
