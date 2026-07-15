# Predict from a `cate_result` object at new targets

Predict from a `cate_result` object at new targets

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

If `interval = "none"` and `se.fit = FALSE`, a numeric vector of
predictions. Otherwise a matrix/list analogous to
[`predict.lm()`](https://rdrr.io/r/stats/predict.lm.html).

## Examples

``` r
set.seed(3)
n <- 40
x_mat <- matrix(rnorm(n * 2), n, 2)
fit <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
predict(fit, x0 = c(0.2, -0.1))
#> [1] 0.191884
```
