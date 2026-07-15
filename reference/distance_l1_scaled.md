# Scaled L1 distance spec (one-time setup + eval/batch)

Builder that computes per-column scales once (`setup(X)`) and reuses
them in [`eval()`](https://rdrr.io/r/base/eval.html) / `eval_batch()`.
Uses a compiled fast path when available.

## Usage

``` r
distance_l1_scaled(
  rule = c("sd", "mad", "iqr", "range", "custom"),
  scale_vec = NULL,
  eps = 1e-12
)
```

## Arguments

- rule:

  Scaling rule: `"sd"`, `"mad"`, `"iqr"`, `"range"`, or `"custom"`.

- scale_vec:

  Optional numeric vector of length \\p\\ for `rule = "custom"`.

- eps:

  Small positive number to avoid division by zero in scaling.

## Value

A list with `name`, `setup()`,
[`eval()`](https://rdrr.io/r/base/eval.html), and `eval_batch()`.

## Examples

``` r
X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l1_scaled(rule = "sd")
ctx <- dspec$setup(X)
dspec$eval_batch(X, X0, ctx)
#>          [,1]     [,2]
#> [1,] 3.219608 2.296217
#> [2,] 5.780489 4.755213
#> [3,] 3.715314 4.143096
#> [4,] 4.970775 5.921876
```
