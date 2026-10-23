# DNN weights (choose method), distance-based ranking

Use `method = "exact"` for the exact finite-population rank-weight rule
and `method = "asymp"` for the asymptotic full-rank approximation.

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

## Details

Fixed top-`k` kernels are available in both cases. The package treats
them as implementation-supported generalizations; the accompanying
manuscript's current derivations are written for `k = 1`.

## Numerical limits

Asymptotic intervals use host Beta probabilities. At subnormal scales,
endpoint rounding can erase an interval or change its zero/nonzero
support; R and the Stata companion can disagree. Signed TDNN
combinations can amplify component rounding. These extreme-tail
limitations remain unresolved. Normalization or switching the R/compiled
deletion engine does not recover probability mass already lost during
table construction. Exact and asymptotic weights are different
constructions; changing the method changes the localization rule.

## Exact distance ties

For each target, observations with exactly equal computed distance share
the average of their assigned rank weights, including distinct covariate
vectors. Equality is exact represented equality, without an epsilon or
rounding. Groups may differ by target. Block mass is preserved for
signed and unnormalized rules. Equal distances are ordered by original
row index. A bare supplied rank or order matrix defines a strict order
and is not averaged. A list containing ranks or orders together with
`dx` preserves that authoritative permutation and averages by the
supplied distances. No outcomes, treatments or nuisance values define
membership.

For a block at ranks `a:b`, each member receives
`sum(w_n[a:b]) / (b-a+1)`. Estimator jackknife deletions remove the
actual observation, update ranks, and average remaining blocks with
`w_(n-1)`. If a member is deleted, each survivor receives
`sum(w_(n-1)[a:(b-1)]) / (b-a)`; a deleted singleton contributes zero.
Exact support edges describe the unaveraged rank table; averaging can
extend nonzero observation weights through a block crossing that edge.
Lower-level rank-only APIs retain their rank-primitive semantics.

A warning is emitted once per public call with distance ties. The
accompanying theory assumes a continuous distance distribution and does
not establish inference under persistent ties, including for `k = 1`.

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0 <- c(0, 0)
out <- weights_dnn(x_mat, x0 = x0, k = 1L, s = 3L, method = "exact")
dim(out$W)
#> [1] 10  1
```
