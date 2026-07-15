# DNN weights (exact full-U), distance-based ranking

Exact DNN weights are derived from the full rank permutation for each
target. They are not constructed by keeping only the first `s`
observations. For the exact finite-population rule, the support in rank
order is `j <= n - s + k`; in particular, for the paper-faithful `k = 1`
DNN this becomes `j <= n - s + 1`. Fixed top-`k` kernels are implemented
for any `1 <= k <= s`. The fixed-`k` extension of the theory is
essentially direct, but the currently written manuscript derivations
cover the `k = 1` case.

## Usage

``` r
weights_dnn_exact(
  X,
  X0 = NULL,
  x0 = NULL,
  distance = "l2",
  k,
  s,
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

- ctx:

  optional context

- return:

  character vector of extras to include; any of
  `c("W","ord","ranks","dx")`

## Value

a `weights_result` with `W` and metadata

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
out <- weights_dnn_exact(x_mat, x0_mat, distance = "l2", k = 1L, s = 3L, return = "W")
dim(out$W)
#> [1] 10  3
```
