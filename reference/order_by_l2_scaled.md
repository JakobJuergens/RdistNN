# Order rows by scaled L2 (stable ties)

Computes per-column scales (or uses `scale_vec` for `"custom"`), then
prefers compiled `rustat_order_by_l2_scaled()` if available. Otherwise
computes distances via
[`dist_l2_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/dist_l2_scaled.md)
and orders deterministically.

## Usage

``` r
order_by_l2_scaled(
  X,
  x0,
  rule = c("sd", "mad", "iqr", "range", "custom"),
  scale_vec = NULL,
  eps = 1e-12,
  squared = FALSE
)
```

## Arguments

- X:

  Numeric matrix \\n \times p\\.

- x0:

  Numeric vector of length \\p\\.

- rule:

  Scaling rule: `"sd"`, `"mad"`, `"iqr"`, `"range"`, or `"custom"`.

- scale_vec:

  Optional numeric vector of length \\p\\ for `rule = "custom"`.

- eps:

  Small positive number to avoid division by zero in scaling.

- squared:

  Logical; if `TRUE`, return squared distances.

## Value

A list with `dx` and `ord`.

## Examples

``` r
X <- matrix(rnorm(10), 5, 2)
order_by_l2_scaled(X, c(0, 0), rule = "sd")
#> $dx
#> [1] 0.8898048 1.0098536 2.1772303 1.2865692 0.4105007
#> 
#> $ord
#> [1] 5 1 2 4 3
#> 
```
