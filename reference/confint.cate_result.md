# Confidence intervals for `cate_result`

Confidence intervals for `cate_result`

## Usage

``` r
# S3 method for class 'cate_result'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `cate_result` object from
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md).

- parm:

  Optional integer vector selecting target indices.

- level:

  Confidence level.

- ...:

  Not used.

## Value

A matrix with lower and upper confidence bounds.

## Examples

``` r
set.seed(2)
n <- 25
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
confint(out)
#>           conf.low conf.high
#> target_1 -1.285905 0.6939418
```
