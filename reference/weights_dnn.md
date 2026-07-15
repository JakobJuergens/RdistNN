# DNN weights (choose method), distance-based ranking

Use `method = "exact"` for the exact finite-population rank-weight rule
and `method = "asymp"` for the asymptotic full-rank approximation. Fixed
top-`k` kernels are available in both cases. The package treats them as
implementation-supported generalizations; the accompanying manuscript's
current derivations are written for `k = 1`.

## Usage

``` r
weights_dnn(
  X,
  X0 = NULL,
  x0 = NULL,
  distance = "l2",
  k,
  s,
  method = c("exact", "asymp"),
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

- method:

  one of `"exact"` or `"asymp"`

- normalize:

  logical; renormalize columns to sum 1

- ctx:

  optional context

- return:

  character vector of extras to include; any of
  `c("W","ord","ranks","dx")`

## Value

A `weights_result` object.

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0 <- c(0, 0)
out <- weights_dnn(x_mat, x0 = x0, k = 1L, s = 3L, method = "exact")
dim(out$W)
#> [1] 10  1
```
