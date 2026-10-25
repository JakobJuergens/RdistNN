#ifndef RDISTNN_RANK_WEIGHT_TABLES_HPP
#define RDISTNN_RANK_WEIGHT_TABLES_HPP
#include "compensated_arithmetic.hpp"

// Portable order-statistic PMFs. Anchoring at a mode avoids enormous
// binomials. A separate binary exponent prevents subnormal recurrence plateaus.
// Callers validate 1 <= k <= s <= n. Workspace is linear in n.
namespace rdistnn_weights {
// Hosts may supply a calling-thread interrupt check. Portable/native callers
// keep the default null callback and unchanged arithmetic.
using interrupt_check = void (*)();
inline int checked_end(int begin, int end, interrupt_check check) {
  if (!check)
    return end;
  check();
  return begin + std::min(65536, end - begin);
}
inline int checked_begin(int end, int begin, interrupt_check check) {
  if (!check)
    return begin;
  check();
  return end - std::min(65536, end - begin);
}
// Keep the callback outside this accumulator loop: inlining it can make GCC
// spill the compensated state on every addition, rather than at chunk boundaries.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
inline rdistnn_arithmetic::accumulator
sum_chunk(const double *values, int n, rdistnn_arithmetic::accumulator total) {
  for (int i = 0; i < n; ++i)
    total.add(values[i]);
  return total;
}
inline double table_sum(const double *values, int n, interrupt_check check) {
  if (!check)
    return rdistnn_arithmetic::sum(values, n);
  rdistnn_arithmetic::accumulator total;
  for (int begin = 0; begin < n;) {
    const int end = checked_end(begin, n, check);
    total = sum_chunk(values + begin, end - begin, total);
    begin = end;
  }
  return total.value();
}
struct scaled_probability {
  double mantissa = .5;
  long long exponent = 1;
  double step(double ratio) {
    mantissa *= ratio;
    while (mantissa < .5) {
      mantissa *= 2;
      --exponent;
    }
    while (mantissa >= 1) {
      mantissa /= 2;
      ++exponent;
    }
    return exponent < -1074 ? 0.0 : mantissa * std::ldexp(1.0, static_cast<int>(exponent));
  }
};
inline std::vector<double> exact_table(int n, int k, int s, interrupt_check check = nullptr) {
  if (check)
    check();
  std::vector<double> out(n, 0.0);
  if (k == s) {
    const double uniform = 1.0 / n;
    for (int begin = 0; begin < n;) {
      const int end = checked_end(begin, n, check);
      std::fill(out.begin() + begin, out.begin() + end, uniform);
      begin = end;
    }
    return out;
  }
  const int width = n - s + 1;
  std::vector<double> component(width), correction(n, 0.0);
  for (int l = 1; l <= k; ++l) {
    if (check)
      check();
    for (int begin = 0; begin < width;) {
      const int end = checked_end(begin, width, check);
      std::fill(component.begin() + begin, component.begin() + end, 0.0);
      begin = end;
    }
    const int last = n - s + l;
    const int mode = std::max(
        l,
        std::min(last,
                 static_cast<int>(std::floor((static_cast<double>(n) * (l - 1)) / (s - 1))) + 1));
    component[mode - l] = 1;
    scaled_probability right, left;
    for (int begin = mode; begin < last;) {
      const int end = checked_end(begin, last, check);
      for (int rank = begin; rank < end; ++rank)
        component[rank - l + 1] = right.step((static_cast<double>(rank) / (rank - l + 1)) *
                                             (static_cast<double>(n - rank - s + l) / (n - rank)));
      begin = end;
    }
    for (int end = mode; end > l;) {
      const int begin = checked_begin(end, l, check);
      for (int rank = end; rank > begin; --rank)
        component[rank - l - 1] =
            left.step((static_cast<double>(rank - l) / (rank - 1)) *
                      (static_cast<double>(n - rank + 1) / (n - rank - s + l + 1)));
      end = begin;
    }
    const double mass = table_sum(component.data(), width, check);
    for (int begin = l - 1; begin < last;) {
      const int end = checked_end(begin, last, check);
      for (int i = begin; i < end; ++i) {
        rdistnn_arithmetic::accumulator total;
        total.high = out[i];
        total.low = correction[i];
        total.add(component[i - l + 1] / mass);
        out[i] = total.high;
        correction[i] = total.low;
      }
      begin = end;
    }
  }
  for (int begin = 0; begin < n;) {
    const int end = checked_end(begin, n, check);
    for (int i = begin; i < end; ++i)
      out[i] = (out[i] + correction[i]) / k;
    begin = end;
  }
  const double mass = table_sum(out.data(), n, check);
  for (int begin = 0; begin < n;) {
    const int end = checked_end(begin, n, check);
    for (int i = begin; i < end; ++i)
      out[i] /= mass;
    begin = end;
  }
  return out;
}
} // namespace rdistnn_weights
#endif
