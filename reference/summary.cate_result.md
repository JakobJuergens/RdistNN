# Summarize a `cate_result` object

Summarize a `cate_result` object

## Usage

``` r
# S3 method for class 'cate_result'
summary(object, level = 0.95, ...)
```

## Arguments

- object:

  A `cate_result` object from
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md).

- level:

  Confidence level used for interval columns.

- ...:

  Not used.

## Value

An object of class `summary.cate_result` with a coefficient table and
`meta`, including the original fit's signal provenance and creation
version.

## Examples

``` r
set.seed(1)
n <- 30
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n))
summary(out)
#> <summary.cate_result>
#>   signal: supplied
#>   localization: s = 5, s/n = 0.167
#>   targets: 1
#>          estimate std.error statistic p.value conf.low conf.high
#> target_1   0.2394    0.2734    0.8758  0.3811  -0.2964    0.7753
```
