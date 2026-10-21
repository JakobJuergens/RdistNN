# Reporting and prediction for nonparametric regression results

Methods for `np_result`, including the `tdnn_np_result` subclass
returned by
[`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md).
Estimates are conditional means, not CATE coefficients.

## Usage

``` r
# S3 method for class 'np_result'
coef(object, ...)

# S3 method for class 'np_result'
vcov(object, ...)

# S3 method for class 'np_result'
summary(object, level = 0.95, ...)

# S3 method for class 'summary.np_result'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'np_result'
confint(object, parm = NULL, level = 0.95, ...)

# S3 method for class 'np_result'
predict(
  object,
  newdata = NULL,
  X0 = NULL,
  x0 = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  An object from
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md)
  or
  [`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md).

- ...:

  Not used.

- level:

  Confidence level in (0,1).

- x:

  A summary object for printing.

- digits:

  Printing precision.

- parm:

  Optional integer vector selecting target indices.

- newdata:

  Optional matrix/data frame of target covariates. Data-frame fits
  extract the retained covariate column names, allowing extra columns.

- X0, x0:

  Alternative target matrix or single target vector.

- se.fit:

  Logical; return jackknife standard errors with predictions.

- interval:

  One of `"none"` or `"confidence"`.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named vector of
local means; [`vcov()`](https://rdrr.io/r/stats/vcov.html) an `m x m`
covariance matrix in target order;
[`confint()`](https://rdrr.io/r/stats/confint.html) a two-column matrix.
[`summary()`](https://rdrr.io/r/base/summary.html) returns class
`summary.np_result` with `coefficients`, `level`, `has_se` and `meta`.
Prediction returns a vector, a `fit`/`se.fit` list, or a
`fit`/`lwr`/`upr` matrix (inside the list when `se.fit = TRUE`).
Printing returns its input invisibly.

## Details

Fits retain canonical outcomes, training covariates, weight controls and
the resolved distance specification and fitted context independently of
optional `return` extras. Prediction uses this full training state and
fitted geometry, including after serialization; it never fits on
newdata. Custom distance functions must respect their retained context
to freeze fitted preprocessing. A named distance is resolved at fitting
time.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the full joint
covariance across targets from common delete-one observations. Its
diagonal agrees with the fit's squared SEs. The `jackknife = "loo"` fit
and covariance recompute each deletion's geometry using the original
context. Prediction uncertainty always uses frozen training geometry and
rank deletions, as does `jackknife = "fast"`. Shared scaled covariance
arithmetic recovers representable results after intermediate overflow or
underflow. Unrepresentable nonzero variances raise a rescaling error
instead of giving infinite or false-zero SEs.

Intervals are pointwise normal-reference intervals. A zero SE gives
collapsed Wald bounds and, for a nonzero mean, an infinite statistic and
zero p-value; zero divided by zero gives undefined statistics/p-values.
Zero estimated variance does not establish negligible bias or coverage.
Missing fitted uncertainty gives `NA` summary columns and an error from
[`confint()`](https://rdrr.io/r/stats/confint.html)/[`vcov()`](https://rdrr.io/r/stats/vcov.html).
Prediction can request uncertainty even if the original fit used
`se = FALSE`, provided the orders permit deletion.

## Examples

``` r
fit <- np_dnn(X = matrix(1:7, ncol = 1), X0 = matrix(c(2, 5), ncol = 1),
              Y = c(2, 1, 4, 3, 6, 5, 7), s = 3, return = character())
coef(fit)
#> target_1 target_2 
#> 2.114286 4.828571 
summary(fit)
#> <summary.np_result>
#>   localization: s = 3, s/n = 0.429
#>          estimate std.error statistic   p.value conf.low conf.high
#> target_1    2.114    0.9726     2.174 2.971e-02   0.2081     4.021
#> target_2    4.829    1.0100     4.782 1.732e-06   2.8500     6.807
vcov(fit)
#>           target_1  target_2
#> target_1 0.9459184 0.1111224
#> target_2 0.1111224 1.0193878
confint(fit)
#>           conf.low conf.high
#> target_1 0.2080574  4.020514
#> target_2 2.8496990  6.807444
predict(fit, x0 = 4, se.fit = TRUE)
#> $fit
#> [1] 3.685714
#> 
#> $se.fit
#> [1] 0.739829
#> 
```
