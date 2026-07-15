# Scaled L2 distance (per-coordinate scaling)

Computes \\\\(X\_{i,\cdot} - x_0) \odot \mathrm{inv\\s}\\\_2\\ (or
squared).

## Usage

``` r
dist_l2_scaled(X, x0, inv_s, squared = FALSE)
```

## Arguments

- X:

  Numeric matrix \\n \times p\\.

- x0:

  Numeric vector of length \\p\\.

- inv_s:

  Numeric vector of length \\p\\ with elementwise scales.

- squared:

  Logical; if `TRUE`, return squared distances.

## Value

Numeric vector of length \\n\\.

## Examples

``` r
X <- matrix(rnorm(12), 4, 3)
dist_l2_scaled(X, c(0, 0, 0), inv_s = c(1, 2, 0.5))
#> [1] 0.6183306 0.6208253 1.7850259 1.4034221
```
