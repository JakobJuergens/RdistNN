# DNN weights (asymptotic full-U), distance-based ranking

The asymptotic DNN rule is evaluated on the full rank grid and is
therefore typically treated as a full-support approximation rather than
as a hard finite-support cutoff. Numerical underflow can still produce
zero entries.

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

## Numerical limits

Asymptotic intervals use host Beta probabilities. At subnormal scales,
endpoint rounding can erase an interval or change its zero/nonzero
support; R and the Stata companion can disagree. Signed TDNN
combinations can amplify component rounding. These extreme-tail
limitations remain unresolved. Normalization or switching the R/compiled
deletion engine does not recover probability mass already lost during
table construction. Exact and asymptotic weights are different
constructions; changing the method changes the localization rule.

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
x0_mat <- matrix(rnorm(4), 2, 2)
out <- weights_dnn_asymp(x_mat, x0_mat, k = 1L, s = 3L, normalize = TRUE)
dim(out$W)
#> [1] 10  2
```
