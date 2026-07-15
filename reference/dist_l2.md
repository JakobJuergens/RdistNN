# Fast L2 (Euclidean) distance to a single query point

Computes \\d_i = \\X\_{i,\cdot} - x_0\\\_2\\ for each row.

## Usage

``` r
dist_l2(X, x0, squared = FALSE)
```

## Arguments

- X:

  Numeric matrix \\n \times p\\.

- x0:

  Numeric vector of length \\p\\.

- squared:

  Logical; if `TRUE`, return squared distances.

## Value

Numeric vector of length \\n\\.

## Examples

``` r
X <- matrix(rnorm(12), 4, 3)
dist_l2(X, c(0, 0, 0))
#> [1] 0.8735564 0.4251124 1.9407667 1.5316326
```
