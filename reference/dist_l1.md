# Fast L1 (Manhattan) distance to a single query point

Computes \\d_i = \\X\_{i,\cdot} - x_0\\\_1\\ for each row
\\i=1,\dots,n\\. If the compiled symbol
[`rustat_dist_l1()`](https://JakobJuergens.github.io/RdistNN/reference/internal_fastpaths.md)
is registered (via Rcpp), it is used; otherwise a vectorised R fallback
is used.

## Usage

``` r
dist_l1(X, x0)
```

## Arguments

- X:

  Numeric matrix \\n \times p\\. Each row is an observation.

- x0:

  Numeric vector of length \\p\\. The query point.

## Value

Numeric vector of length \\n\\ with L1 distances.

## See also

[`dist_l2()`](https://JakobJuergens.github.io/RdistNN/reference/dist_l2.md),
[`dist_l1_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/dist_l1_scaled.md),
[`distance_l1()`](https://JakobJuergens.github.io/RdistNN/reference/distance_l1.md),
[`order_by_l1()`](https://JakobJuergens.github.io/RdistNN/reference/order_by_l1.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(12), 4, 3)
x0 <- c(0, 1, -1)
dist_l1(X, x0)
#> [1] 2.872727 2.698723 3.859981 3.246799
```
