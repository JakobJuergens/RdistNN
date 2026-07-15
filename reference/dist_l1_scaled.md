# Scaled L1 distance (per-coordinate scaling)

Computes \\\\(X\_{i,\cdot} - x_0) \odot \mathrm{inv\\s}\\\_1\\ for each
row. If `rustat_dist_l1_scaled()` exists it is used; else inputs are
rescaled in R and passed to
[`dist_l1()`](https://JakobJuergens.github.io/RdistNN/reference/dist_l1.md).

## Usage

``` r
dist_l1_scaled(X, x0, inv_s)
```

## Arguments

- X:

  Numeric matrix \\n \times p\\. Each row is an observation.

- x0:

  Numeric vector of length \\p\\. The query point.

- inv_s:

  Numeric vector of length \\p\\ with elementwise scales (often the
  reciprocals of column scale estimates).

## Value

Numeric vector of length \\n\\.

## Examples

``` r
X <- matrix(rnorm(12), 4, 3)
dist_l1_scaled(X, c(0, 0, 0), inv_s = c(1, 2, 0.5))
#> [1] 1.113110 4.493440 2.804656 2.227412
```
