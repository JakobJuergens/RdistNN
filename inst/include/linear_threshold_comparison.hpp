#ifndef RDISTNN_LINEAR_THRESHOLD_COMPARISON_HPP
#define RDISTNN_LINEAR_THRESHOLD_COMPARISON_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace rdistnn_arithmetic {
namespace threshold_detail {

static_assert(std::numeric_limits<double>::radix == 2 &&
                  std::numeric_limits<double>::digits == 53 &&
                  std::numeric_limits<double>::min_exponent == -1021 &&
                  std::numeric_limits<double>::max_exponent == 1024,
              "PLQR comparisons require binary64 arithmetic");

struct binary_operand {
  double mantissa;
  int exponent;
};

inline binary_operand decompose(double value) {
  if (value == 0)
    return {0.0, 0};
  const int power = std::max(-1074, std::min(1022, std::ilogb(value)));
  return {std::ldexp(value, -power), power};
}

// Every limb is an exact nonnegative integer below 2^26. A binary64
// significand has 53 bits; the exponent-1022 cap permits one extra bit,
// which is zero at the low end. Three limbs therefore hold its exact value.
constexpr int digit_bits = 26;
constexpr int fraction_bits = std::numeric_limits<double>::digits - 1;
constexpr int operand_bits = std::numeric_limits<double>::digits + 1;
constexpr int operand_chunks = (operand_bits + digit_bits - 1) / digit_bits;
constexpr double limb_base = 67108864.0; // 2^digit_bits
constexpr int minimum_operand_power =
    std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits;
constexpr int bit_origin = 2 * minimum_operand_power - 2 * fraction_bits;
// Each product is below 2^2048. At most INT_MAX+2 terms make either sign's
// subtotal smaller than 2^2080 on the supported 32-bit matrix dimensions.
// Origin -2252, 170 base-2^26 limbs cover through exponent2167, including
// three extra carry limbs. No product bits or subnormal residuals are dropped.
constexpr int product_exclusive_power = 2 * std::numeric_limits<double>::max_exponent;
constexpr int term_count_power = std::numeric_limits<int>::digits + 1;
constexpr int required_bits = product_exclusive_power + term_count_power - bit_origin;
constexpr std::size_t limb_count = (required_bits + digit_bits - 1) / digit_bits + 3;
using exact_limbs = std::array<double, limb_count>;

inline bool add_piece(exact_limbs &target, double piece, int bit) {
  if (piece == 0)
    return true;
  if (bit < 0)
    return false;
  std::size_t index = static_cast<std::size_t>(bit / digit_bits);
  piece = std::ldexp(piece, bit % digit_bits);
  while (piece != 0) {
    if (index >= target.size())
      return false;
    // Split products before shifting: piece<2^51, so this integer addition
    // is below 2^52 and exact. Following carry additions are below 2^27.
    const double next = target[index] + piece;
    const double carry = std::floor(next / limb_base);
    target[index] = next - carry * limb_base;
    piece = carry;
    ++index;
  }
  return true;
}

inline bool add_exact_product(exact_limbs &positive, exact_limbs &negative, double left,
                              double right) {
  if (!std::isfinite(left) || !std::isfinite(right))
    return false;
  if (left == 0 || right == 0)
    return true;
  const auto a = decompose(left), b = decompose(right);
  double x = std::ldexp(std::abs(a.mantissa), fraction_bits);
  double y = std::ldexp(std::abs(b.mantissa), fraction_bits);
  std::array<double, operand_chunks> xs{}, ys{};
  for (int i = 0; i < operand_chunks; ++i) {
    const double next_x = std::floor(x / limb_base), next_y = std::floor(y / limb_base);
    xs[i] = x - next_x * limb_base;
    ys[i] = y - next_y * limb_base;
    x = next_x;
    y = next_y;
  }
  if (x != 0 || y != 0)
    return false;
  auto &target = ((left < 0) == (right < 0)) ? positive : negative;
  const int start = a.exponent + b.exponent - 2 * fraction_bits - bit_origin;
  for (int i = 0; i < operand_chunks; ++i) {
    for (int j = 0; j < operand_chunks; ++j) {
      // Each product of two 26-bit integers is below 2^52 and exact.
      const double product = xs[i] * ys[j];
      const double high = std::floor(product / limb_base), low = product - high * limb_base;
      const int bit = start + digit_bits * (i + j);
      if (!add_piece(target, low, bit) || !add_piece(target, high, bit + digit_bits))
        return false;
    }
  }
  return true;
}

// Return the exact sign -1/0/+1, or 2 for invalid internal operands/state.
inline int exact_sign(const double *x, std::ptrdiff_t stride, const double *coefficients, int p,
                      double g, double y) {
  exact_limbs positive{}, negative{};
  for (int j = 0; j < p; ++j)
    if (!add_exact_product(positive, negative, x[static_cast<std::ptrdiff_t>(j) * stride],
                           coefficients[j]))
      return 2;
  // Keeping y*(-1) as a product also avoids forming an out-of-range -y
  // in hosts whose positive finite range is smaller than the negative range.
  if (!add_exact_product(positive, negative, g, 1.0) ||
      !add_exact_product(positive, negative, y, -1.0))
    return 2;
  for (std::size_t i = limb_count; i-- > 0;) {
    if (positive[i] > negative[i])
      return 1;
    if (positive[i] < negative[i])
      return -1;
  }
  return 0;
}

} // namespace threshold_detail

// Compare y <= dot(x,coefficients)+g. Preserve a decisively separated ordinary
// comparison; otherwise determine its exact sign from the original binary64
// operands. Return -1 only for invalid internal operands/state, otherwise0/1.
// The helper owns bounded scalar workspace and calls no host API.
inline int linear_threshold_indicator(const double *x, std::ptrdiff_t stride,
                                      const double *coefficients, int p, double g, double y,
                                      double ordinary_threshold) {
  if (p < 0 || (p > 0 && (!x || !coefficients)))
    return -1;
  double mass = std::abs(g) + std::abs(y);
  for (int j = 0; j < p; ++j) {
    volatile double product = x[static_cast<std::ptrdiff_t>(j) * stride] * coefficients[j];
    mass += std::abs(product);
  }
  const double margin = ordinary_threshold - y;
  if (std::isfinite(ordinary_threshold) && std::isfinite(margin) && std::isfinite(mass) &&
      mass < std::ldexp(1.0, 1022)) {
    const double count = static_cast<double>(p) + 2.0;
    // With N=2*count rounded operations and u=eps/2, gamma_N<2*N*u.
    // Accounting also for rounding the nonnegative mass gives error less
    // than 4*count*eps*mass+5*count*eta. The factors below include slack for
    // the filter's own rounding. Assumptions are IEEE nearest rounding and
    // gradual underflow; ordinary FMA/wider accumulation is no less accurate.
    volatile double relative = 8.0 * (count + 1.0) * std::numeric_limits<double>::epsilon() * mass;
    volatile double underflow = 32.0 * count * std::numeric_limits<double>::denorm_min();
    const double uncertainty = relative + underflow;
    if (std::abs(margin) > uncertainty)
      return y <= ordinary_threshold ? 1 : 0;
  }
  const int sign = threshold_detail::exact_sign(x, stride, coefficients, p, g, y);
  return sign == 2 ? -1 : (sign >= 0 ? 1 : 0);
}

} // namespace rdistnn_arithmetic
#endif
