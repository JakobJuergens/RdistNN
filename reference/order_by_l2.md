# Order rows by L2 distance to a point (stable ties)

Order rows by L2 distance to a point (stable ties)

## Usage

``` r
order_by_l2(X, x0, squared = FALSE)
```

## Arguments

- X:

  Numeric matrix \\n \times p\\.

- x0:

  Numeric vector of length \\p\\.

- squared:

  Logical; if `TRUE`, return squared distances.

## Value

A list with `dx` and `ord`.

## Examples

``` r
X <- matrix(rnorm(10), 5, 2)
order_by_l2(X, c(0, 0))
#> $dx
#> [1] 1.7287137 2.0749119 1.2744925 0.4685352 2.2873867
#> 
#> $ord
#> [1] 4 3 1 2 5
#> 
```
