# Nonparametric regression **singleton** score from observed outcome \\Y\\ (order-1 pure score; no weights)

Constructs a **pure-score** kernel for nonparametric regression that
uses only the observed response \\Y\\. For a singleton subset
\\S=\\i\\\\ the kernel returns exactly \\h(S) = Y_i\\. This score is
**order-1 only** (s = 1). It performs no nearest-neighbor ordering,
kernel weighting, or L-stat aggregation. Any weighting/averaging of
subset evaluations is handled outside the score layer (e.g., by
rank-/distance-weight modules).

## Usage

``` r
score_np_y(x0, x_cols, s = 1, y_col = NULL, y = NULL)
```

## Arguments

- x0:

  Numeric length-\\p\\ target point. Used only to check that
  `length(x0)` matches `length(x_cols)`; it is not used in the kernel.

- x_cols:

  Integer or character indices of the \\X\\ columns in `data`, used for
  dimension checks with `x0`.

- s:

  Intended subset size (kernel order). **Must be 1**; otherwise an error
  is thrown.

- y_col:

  Optional name of the column in `data` that contains the response
  \\Y\\. If `NULL`, supply `y`.

- y:

  Optional numeric vector (length \\n\\) of \\Y\\ values to use when
  `y_col` is `NULL`.

## Value

An `ustat_score` object with:

- `$name`: a descriptive name,

- `$order`: always `1L`,

- `$setup`: stores `g <- Y` in the context,

- `$eval`/`$eval_batch`: return \\Y_i\\ for each singleton subset.

## Details

The arguments `x0`, `x_cols`, and `s` are accepted for API consistency
and basic dimension checks; only \\Y\\ is used to construct the score.
This factory **enforces** `s = 1`.

## Examples

``` r
set.seed(1)
n <- 50
p <- 2
dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
dat$Y <- with(dat, X1 - 2 * X2 + rnorm(n))

sc <- score_np_y(
  x0 = c(0, 0),
  x_cols = c("X1", "X2"),
  s = 1L,
  y_col = "Y"
)
ctx <- score_setup(sc, dat)
subs <- matrix(seq_len(n), ncol = 1L) # singleton subsets
score_eval_batch(sc, subs, dat, ctx) # returns dat$Y
#>  [1] -2.043032249  1.449811984 -2.428789644  4.012035767 -3.191124276
#>  [6] -3.013980912  1.938579481  3.736768187 -0.179472545  1.646896901
#> [11] -3.927190807  0.006678512 -0.568437067 -2.921400558  2.404096592
#> [16] -0.815326137  3.273734126 -2.266386815  1.008902850 -3.928652502
#> [21] -0.537999148  3.545067988 -1.361467132 -0.300712963  3.026901807
#> [26]  0.073645096  0.657223836 -1.510597259 -1.308493182  1.272713180
#> [31]  2.556177458 -0.421324983 -1.437006189  1.474934478 -2.258394071
#> [36] -2.617345129 -2.821465755  0.020774546 -0.292107029  0.172081389
#> [41] -0.993842960 -1.492513980 -3.288814292 -1.307294502 -4.978342709
#> [46] -2.575287009  5.004932925  1.932459373  1.050578487  0.187303465
```
