# Scaled L2 distance spec (one-time setup + eval/batch)

Builder that computes per-column scales once (`setup(X)`) and reuses
them in [`eval()`](https://rdrr.io/r/base/eval.html) / `eval_batch()`.
Uses a compiled fast path when available.

## Usage

``` r
distance_l2_scaled(
  rule = c("sd", "mad", "iqr", "range", "custom"),
  scale_vec = NULL,
  eps = 1e-12,
  squared = FALSE
)
```

## Arguments

- rule:

  Scaling rule: `"sd"`, `"mad"`, `"iqr"`, `"range"`, or `"custom"`.

- scale_vec:

  Optional numeric vector of length \\p\\ for `rule = "custom"`.

- eps:

  Small positive number to avoid division by zero in scaling.

- squared:

  Logical; if `TRUE`, return squared distances.

## Value

A list with `name`, `setup()`,
[`eval()`](https://rdrr.io/r/base/eval.html), and `eval_batch()`.

## Examples

``` r
X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l2_scaled(rule = "sd")
ctx <- dspec$setup(X)
dspec$eval_batch(X, X0, ctx)
#>           [,1]     [,2]
#> [1,] 2.8710813 2.148632
#> [2,] 0.8411369 1.666764
#> [3,] 2.9497331 1.985105
#> [4,] 3.1551688 1.698267
```
