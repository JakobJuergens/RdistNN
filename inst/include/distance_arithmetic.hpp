#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace rustat {

// Stata reserves this range for missing values. Share exceptional subtraction
// selection with its native/Mata engines even when IEEE subtraction is finite.
inline constexpr double distance_range_limit = 0x1p1023;

inline double scaled_distance_difference(double x, double target, double inverse) {
  const double difference = x - target;
  if (std::fabs(difference) >= distance_range_limit && std::isfinite(x) && std::isfinite(target) &&
      std::fabs(inverse) < 1.0) {
    const volatile double half_x = x * 0.5, half_target = target * 0.5;
    const double half_difference = half_x - half_target;
    if (std::fabs(half_difference) >= distance_range_limit) {
      const volatile double quarter_x = x * 0.25, quarter_target = target * 0.25;
      return (quarter_x - quarter_target) * (4.0 * inverse);
    }
    return half_difference * (inverse + inverse);
  }
  return difference * inverse;
}

// A split power keeps a representable final subnormal from rounding twice.
inline double distance_power2(double value, int power) {
  if (value == 0.0)
    return value;
  const int half = power / 2;
  const volatile double first = value * std::ldexp(1.0, half);
  return first * std::ldexp(1.0, power - half);
}

// Exceptional scaled rows retain coordinate products outside the double
// exponent range until the complete norm is restored. Ordinary rows never
// enter this computation. Input and target arrays are column-major.
inline double scaled_norm_row(const double *X, int n, int p, int i, const double *target,
                              int stride_target, const double *inverse, bool l2,
                              bool squared = false) {
  double scale_m = 0.0, total = 0.0;
  double rounded_max = 0.0;
  int scale_e = 0, active = 0;
  for (int j = 0; j < p; ++j) {
    const double x = X[i + static_cast<std::size_t>(j) * n];
    const double t = target[static_cast<std::size_t>(j) * stride_target];
    double raw = x - t;
    int shift = 0;
    if (std::fabs(raw) >= distance_range_limit) {
      const volatile double quarter_x = x * 0.25, quarter_t = t * 0.25;
      raw = quarter_x - quarter_t;
      shift = 2;
    }
    if (raw == 0.0 || inverse[j] == 0.0)
      continue;
    if (!std::isfinite(raw) || !std::isfinite(inverse[j]))
      return std::fabs(raw * inverse[j]);
    ++active;
    rounded_max = std::max(rounded_max, std::fabs(scaled_distance_difference(x, t, inverse[j])));
    int raw_e, inverse_e;
    const double raw_m = std::frexp(std::fabs(raw), &raw_e);
    const double inverse_m = std::frexp(std::fabs(inverse[j]), &inverse_e);
    volatile double product = raw_m * inverse_m;
    double m = product;
    int e = raw_e + inverse_e + shift;
    if (m < 0.5) {
      m *= 2.0;
      --e;
    }
    if (scale_m == 0.0) {
      scale_m = m;
      scale_e = e;
      total = 1.0;
    } else if (e > scale_e || (e == scale_e && m > scale_m)) {
      const double ratio = distance_power2(scale_m / m, scale_e - e);
      const volatile double square = ratio * ratio;
      const volatile double term = total * (l2 ? square : ratio);
      total = 1.0 + term;
      scale_m = m;
      scale_e = e;
    } else {
      const double ratio = distance_power2(m / scale_m, e - scale_e);
      const volatile double square = ratio * ratio;
      total += l2 ? square : ratio;
    }
  }
  if (scale_m == 0.0)
    return 0.0;
  // A one-coordinate norm is its directly rounded product. Reconstructing
  // it through a normal mantissa can double-round at subnormal boundaries.
  if (active == 1) {
    const volatile double square = rounded_max * rounded_max;
    return squared ? square : rounded_max;
  }
  if (squared) {
    const volatile double first = scale_m * total;
    return distance_power2(first * scale_m, 2 * scale_e);
  }
  const double norm = distance_power2(scale_m * (l2 ? std::sqrt(total) : total), scale_e);
  // Monotonic rounding preserves the exact norm's lower bound at each
  // coordinate, even when normalized-product rounding loses a boundary bit.
  return std::max(norm, rounded_max);
}

} // namespace rustat
