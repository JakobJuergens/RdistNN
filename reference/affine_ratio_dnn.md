# DNN inference for a ratio of two supplied signal averages

Localizes two aligned, externally supplied signals with common DNN
weights. This estimator uses paired frozen-score covariance; it never
divides individual signals or deleted component averages.

## Usage

``` r
affine_ratio_dnn(
  data = NULL,
  x_cols = NULL,
  numerator_col = NULL,
  denominator_col = NULL,
  X = NULL,
  X0 = NULL,
  x0 = NULL,
  numerator = NULL,
  denominator = NULL,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L,
  s = NULL,
  tuning = NULL,
  tuning_mode = c("default", "grid"),
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  denominator_tol = sqrt(.Machine$double.eps),
  row_ids = NULL,
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- data:

  Optional data frame or tibble. Do not combine with X or direct sample
  vectors.

- x_cols:

  Reporting-covariate columns in data, in target-coordinate order.

- numerator_col, denominator_col:

  Signal columns for the data-frame route.

- X:

  Numeric reporting matrix, observations in rows.

- X0:

  Numeric target matrix, targets in rows.

- x0:

  Numeric vector for a single target. Supply exactly one of X0,x0.

- numerator, denominator:

  Finite numeric signal vectors for the matrix route.

- distance:

  Existing distance specification, registered name or function. Geometry
  is frozen across deletions. Exact distance ties use original
  covariates lexicographically; identical covariate rows use input
  position.

- method:

  One of "exact" or "asymp", as in
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
  and
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md).

- k:

  Number of neighbors averaged within each size-s subsample, 1 \<= k \<=
  s.

- s:

  Deterministic integer order. Covariance requires s \< n.

- tuning:

  Optional existing tuning specification, mutually exclusive with s. If
  neither is provided, the package paper rule resolves the order.

- tuning_mode:

  "default" retains a single fit; "grid" also evaluates feasible tuning
  candidates using the same frozen signals and ranks.

- normalize:

  Logical normalization control passed to the asymptotic weight rule, as
  in
  [`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md).
  Exact weights are normalized by construction.

- se:

  Whether to compute covariance.

- jackknife:

  "fast" uses the existing prefix/suffix engine; "loo" uses explicit
  frozen-rank deletions; "none" disables covariance.

- denominator_tol:

  Nonnegative numerical cancellation tolerance. Changes a diagnostic
  flag only, never the ratio or covariance.

- row_ids:

  Optional aligned unique character/numeric identities.

- ctx:

  Optional context passed to the existing distance machinery.

- return:

  Optional extras: W, ranks, ord, dx, ess, signals, deletions,
  pseudo_values, neighbors or tuning. Unknown values warn and are
  ignored.

## Value

An affine_ratio_result with theta, se, m-by-2 components, lists
component_vcov, theta_vcov and theta_omega (NULL without covariance),
x0, neff, per-target status, row_ids, meta and retained frozen model
state. Optional deletions/pseudo_values have dimensions n-by-2-by-m.
Covariance matrices are marginal across targets, not cross-target
covariances.

## Details

Both full and deleted averages use the same order s. The target and
nuisance predictions remain frozen. If V is the paired estimator-scale
component covariance, the ratio variance is \$\$(V\_{YY}+\hat\theta^2
V\_{DD}-2\hat\theta V\_{YD})/\hat\Delta_D^2.\$\$ The computation uses
centered deletion scores to retain cross-covariance. theta_vcov is
already on the estimator scale; theta_omega multiplies it by neff. Do
not divide theta_vcov by neff again.

Signed signals and finite negative denominators are retained. At exact
zero denominator the components remain available and the ratio is NA.
Computational statuses do not verify causal assumptions, learner
stability, point-target undersmoothing or weak-instrument-uniform
validity. If nonzero component variation underflows its covariance,
uncertainty is unavailable with numerical_failure status; it is not
treated as deterministic.

Prediction retains the resolved distance specification and fitted
context; re-registering the distance name does not change an existing
fit.

## Examples

``` r
fit <- affine_ratio_dnn(
  X = matrix(1:4, ncol = 1), x0 = 0,
  numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
)
fit$theta
#> [1] 4
fit$theta_vcov
#> [[1]]
#>       theta
#> theta    34
#> 
d <- data.frame(x = 1:4, y = c(2, 1, 4, 3), d = c(1, -1, 2, 1))
affine_ratio_dnn(
  data = d, x_cols = "x", numerator_col = "y",
  denominator_col = "d", x0 = 0, s = 2, se = FALSE
)$theta
#> [1] 4
```
