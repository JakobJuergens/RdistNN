# Print a `rdistnn_weights_result` object

Print a `rdistnn_weights_result` object

## Usage

``` r
# S3 method for class 'rdistnn_weights_result'
print(x, ...)
```

## Arguments

- x:

  A `rdistnn_weights_result` object created by
  [`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
out <- weights_dnn_exact(x_mat, x0 = c(0, 0), k = 1L, s = 1L)
print(out)
#> <rdistnn_weights_result>
#>   data  : n=10, p=2, m=1
#>   weight: rank (dnn_exact[k=1,s=1])
```
