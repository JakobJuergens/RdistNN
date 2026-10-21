# Predict from a `cate_result` object at new targets

Uses retained observations, supplied signals and fitted geometry,
including after serialization and without optional `return` extras.
Prediction uncertainty uses frozen-geometry rank deletions, even if the
fit used `jackknife = "loo"`. It can be requested from a fit with
`se = FALSE` when `s < n`. Prediction localizes the retained signal
without establishing new inferential guarantees.

## Usage

``` r
# S3 method for class 'cate_result'
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

  A `cate_result` object from
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md).

- newdata:

  Optional matrix/data frame of target covariates.

- X0:

  Optional numeric matrix of targets (rows).

- x0:

  Optional numeric vector for a single target.

- se.fit:

  Logical; if `TRUE`, also return jackknife standard errors.

- interval:

  One of `"none"` or `"confidence"`.

- level:

  Confidence level for interval output.

- ...:

  Not used.

## Value

A numeric vector when no uncertainty is requested; a list with `fit` and
`se.fit` when `se.fit = TRUE`. With `interval = "confidence"`, `fit` is
a matrix with columns `fit`, `lwr`, `upr`, returned directly unless
`se.fit = TRUE`.

## Examples

``` r
set.seed(3)
n <- 40
x_mat <- matrix(rnorm(n * 2), n, 2)
fit <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
predict(fit, x0 = c(0.2, -0.1))
#> [1] 0.191884
```
