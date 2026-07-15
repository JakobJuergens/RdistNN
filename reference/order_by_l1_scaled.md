# Order rows by scaled L1 (stable ties)

Computes per-column scales (or uses `scale_vec` for `"custom"`), then
prefers compiled `rustat_order_by_l1_scaled()` if available. Otherwise
computes distances via
[`dist_l1_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/dist_l1_scaled.md)
and orders with a deterministic tie rule.

## Usage

``` r
order_by_l1_scaled(
  X,
  x0,
  rule = c("sd", "mad", "iqr", "range", "custom"),
  scale_vec = NULL,
  eps = 1e-12
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

## Value

A list with:

- `dx`: numeric vector (length \\n\\)

- `ord`: integer vector of 1-based row indices

## See also

[`order_by_l1()`](https://JakobJuergens.github.io/RdistNN/reference/order_by_l1.md),
[`dist_l1_scaled()`](https://JakobJuergens.github.io/RdistNN/reference/dist_l1_scaled.md)

## Examples

``` r
set.seed(6)
X <- matrix(rnorm(10), 5, 2)
x0 <- c(0, 0)
order_by_l1_scaled(X, x0, rule = "sd")
#> $dx
#> [1] 0.7131843 2.1683489 1.7983424 1.9852135 1.1982833
#> 
#> $ord
#> [1] 1 5 3 4 2
#> 
```
