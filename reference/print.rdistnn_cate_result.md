# Print a `rdistnn_cate_result` object

Print a `rdistnn_cate_result` object

## Usage

``` r
# S3 method for class 'rdistnn_cate_result'
print(x, ...)
```

## Arguments

- x:

  An object returned by
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
n <- 20
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(
  X = x_mat,
  x0 = c(0, 0),
  pseudo_outcome = rnorm(n),
  se = FALSE
)
print(out)
#> <rdistnn_cate_result>
#>   signal: supplied
#>   localization: s = 3, s/n = 0.15
#>   data  : n=20, p=2, m=1
#>   DNN   : method=exact, k=1, distance=l2
#>   est   : tau[1:min(5,m)] = 0.1451  (SEs not requested)
```
