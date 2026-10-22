# DNN weights (exact full-U), distance-based ranking

Exact DNN weights are derived from the full rank permutation for each
target. They are not constructed by keeping only the first `s`
observations. For the exact finite-population rule, the unaveraged rank
table has nominal support `j <= n - s + k`; for `k = 1` this becomes
`j <= n - s + 1`. Duplicate-block averaging can extend nonzero
observation weights beyond this nominal edge. Fixed top-`k` kernels are
implemented for `1 <= k <= s`. The accompanying
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
gives the theory for the unmodified `k = 1` rule.

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

## Duplicate-covariate blocks

Rows identical on every supplied raw covariate receive the average of
their assigned rank weights, preserving each block's total mass. This
includes caller-supplied ranks and signed or unnormalized rank rules.
Distinct covariate vectors at equal distance remain separate. A warning
identifies this permutation-invariant modification as outside the
paper's theory, including for `k = 1`. Estimator jackknife deletions
repeat the average over each remaining duplicate block using the
selected rule's sample-size-n-1 table. For a block at consecutive ranks
`a:b`, the full weight is `sum(w_n[a:b]) / (b-a+1)` per member; after
deleting a member, each survivor receives
`sum(w_(n-1)[a:(b-1)]) / (b-a)`. A deleted singleton contributes zero.
The smaller block size must be used in the deletion. An exact support
edge describes the unaveraged rank table; a straddling duplicate block
can extend nonzero observation weights beyond it. Lower-level rank-only
APIs retain their rank-primitive semantics.

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
out <- weights_dnn_exact(x_mat, x0_mat, distance = "l2", k = 1L, s = 3L, return = "W")
dim(out$W)
#> [1] 10  3
```
