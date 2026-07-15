# DNN weights (asymptotic full-U), distance-based ranking

The asymptotic DNN rule is evaluated on the full rank grid and is
therefore typically treated as a full-support approximation rather than
as a hard finite-support cutoff.

## Usage

``` r
weights_dnn_asymp(
  X,
  X0 = NULL,
  x0 = NULL,
  distance = "l2",
  k,
  s,
  normalize = TRUE,
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- X:

  numeric matrix (n x p)

- X0:

  numeric matrix (m x p); mutually exclusive with `x0`

- x0:

  numeric length-p single target (convenience for m=1)

- distance:

  distance spec or registered name (e.g., "l2")

- k:

  integer; neighbors used in each s-subset (`1 <= k <= s`)

- s:

  integer; U-statistic order / subset size

- normalize:

  logical; renormalize columns to sum 1

- ctx:

  optional context

- return:

  character vector of extras to include; any of
  `c("W","ord","ranks","dx")`

## Value

A `weights_result` object (same structure as
[`weights_dnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn_exact.md)).

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(4), 2, 2)
out <- weights_dnn_asymp(x_mat, x0_mat, k = 1L, s = 3L, normalize = TRUE)
dim(out$W)
#> [1] 10  2
```
