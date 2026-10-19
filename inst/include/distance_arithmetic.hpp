#pragma once
#include <cmath>

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
    return half_difference * (inverse + inverse);
  }
  return difference * inverse;
}

} // namespace rustat
