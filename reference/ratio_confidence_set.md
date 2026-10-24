# Confidence sets for a ratio of paired component estimates

Inverts the pointwise normal-reference moment inequality without
dividing by an estimated denominator. This experimental helper accepts
the components and paired estimator-scale covariance returned by
[`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md).

## Usage

``` r
ratio_confidence_set(
  components,
  component_vcov,
  level = 0.95,
  tol = 64 * .Machine$double.eps
)

# S3 method for class 'rdistnn_ratio_confidence_sets'
print(x, digits = max(3L, getOption("digits") - 2L), ...)
```

## Arguments

- components:

  Finite numeric m-by-2 matrix, with numerator in the first column and
  denominator in the second; m must be positive.

- component_vcov:

  List of m finite symmetric positive semidefinite 2-by-2 covariance
  matrices on the estimator scale. Singular matrices are allowed.

- level:

  Confidence level strictly between zero and one.

- tol:

  Finite nonnegative relative numerical tolerance. It controls
  covariance symmetry/PSD checks and near-zero discriminant diagnostics.

- x:

  A rdistnn_ratio_confidence_sets object.

- digits:

  Printing precision.

- ...:

  Not used.

## Value

A rdistnn_ratio_confidence_sets list with one record per target. Records
contain target_id, type, intervals (r-by-2 lower/upper matrix),
quadratic, discriminant, numerical_status, level and numerical scaling
metadata. Printing returns its input invisibly.

## Details

For candidate theta, invert \$\$(\hat\Delta_Y-\theta\hat\Delta_D)^2 \le
z^2(V\_{YY}-2\theta V\_{YD}+\theta^2V\_{DD}).\$\$ Inputs must retain
numerator-denominator pairing. Do not divide the supplied covariance by
effective sample size again.

Each set can be bounded, a singleton, two rays, a half-line, all real
numbers or empty. Finite endpoints are included; infinities are open in
printed interval notation. An unresolved numerical result has no
intervals and is distinguished from a mathematically empty set by its
type/status. Near-boundary calculations are reported rather than
rounding coefficients or discriminants to zero. Materially indefinite
covariance is rejected; a negative eigenvalue within tolerance gives an
unresolved result. Deterministic-denominator intervals are computed
directly from the components to preserve small positive uncertainty.
Repeated roots of a rounded polynomial alone do not certify zero
uncertainty; such cases are unresolved unless a component-level
zero-variance identity establishes it.

Power-of-two coordinate scaling avoids intermediate overflow. The
numerical record retains original quadratic coefficients where
representable, plus scaled coefficients for t = theta \*
2^(-ratio_scale_exponent). Multiplying that polynomial by
2^quadratic_scale_exponent gives the polynomial of the coordinate-scaled
components. This metadata permits reconstruction when the original
coefficients exceed floating-point range.

At a fixed true candidate, inclusion follows from validity of that
candidate's local score test under its design, learner, bias and
positive variance conditions. This calculation establishes no uniform
coverage over increasingly weak instruments. Algebraic zero-variance
branches do not automatically satisfy a positive-variance test theorem.

## Examples

``` r
f <- affine_ratio_dnn(
  X = matrix(1:4, ncol = 1), x0 = 0,
  numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
)
ratio_confidence_set(f$components, f$component_vcov)
#> Paired ratio confidence sets
#>   target 1: two_rays (-Inf, -0.60543] U [1.3053, Inf) (ok)
# An exactly zero denominator can give two unbounded rays.
ratio_confidence_set(matrix(c(3, 0), 1), list(diag(2)))
#> Paired ratio confidence sets
#>   target 1: two_rays (-Inf, -1.1588] U [1.1588, Inf) (ok)
```
