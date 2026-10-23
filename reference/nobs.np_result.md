# Number of training observations

Number of training observations

## Usage

``` r
# S3 method for class 'np_result'
nobs(object, ...)

# S3 method for class 'cate_result'
nobs(object, ...)
```

## Arguments

- object:

  A fitted result from
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md)
  or
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md).

- ...:

  Not used.

## Value

The integer number of training observations. This is distinct from the
number of reporting profiles and the local effective sample size.

## Examples

``` r
fit <- np_dnn(X = matrix(1:12), Y = sin(1:12), x0 = 4, s = 3)
#> Warning: Exact distance ties detected; weights are averaged within distance blocks, including jackknife deletions. The accompanying theory assumes a continuous distance distribution and does not establish inference under persistent ties.
nobs(fit)
#> [1] 12
```
