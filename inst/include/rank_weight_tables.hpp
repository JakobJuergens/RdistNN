#ifndef RDISTNN_RANK_WEIGHT_TABLES_HPP
#define RDISTNN_RANK_WEIGHT_TABLES_HPP
#include "compensated_arithmetic.hpp"

// Portable order-statistic PMFs. Anchoring at a mode avoids enormous
// binomials. A separate binary exponent prevents subnormal recurrence plateaus.
// Callers validate 1 <= k <= s <= n. Workspace is linear in n.
namespace rdistnn_weights {
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
inline std::vector<double> exact_table(int n, int k, int s) {
  std::vector<double> out(n, 0.0);
  if (k == s) {
    std::fill(out.begin(), out.end(), 1.0 / n);
    return out;
  }
  const int width = n - s + 1;
  std::vector<double> component(width), correction(n, 0.0);
  for (int l = 1; l <= k; ++l) {
    std::fill(component.begin(), component.end(), 0.0);
    const int last = n - s + l;
    const int mode = std::max(
        l,
        std::min(last,
                 static_cast<int>(std::floor((static_cast<double>(n) * (l - 1)) / (s - 1))) + 1));
    component[mode - l] = 1;
    scaled_probability right, left;
    for (int rank = mode; rank < last; ++rank)
      component[rank - l + 1] = right.step((static_cast<double>(rank) / (rank - l + 1)) *
                                           (static_cast<double>(n - rank - s + l) / (n - rank)));
    for (int rank = mode; rank > l; --rank)
      component[rank - l - 1] =
          left.step((static_cast<double>(rank - l) / (rank - 1)) *
                    (static_cast<double>(n - rank + 1) / (n - rank - s + l + 1)));
    const double mass = rdistnn_arithmetic::sum(component.data(), width);
    for (int i = l - 1; i < last; ++i) {
      rdistnn_arithmetic::accumulator total;
      total.high = out[i];
      total.low = correction[i];
      total.add(component[i - l + 1] / mass);
      out[i] = total.high;
      correction[i] = total.low;
    }
  }
  for (int i = 0; i < n; ++i)
    out[i] = (out[i] + correction[i]) / k;
  const double mass = rdistnn_arithmetic::sum(out.data(), n);
  for (double &value : out)
    value /= mass;
  return out;
}
} // namespace rdistnn_weights
#endif
