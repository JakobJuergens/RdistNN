# TDNN weights (choose method, k = 1), distance-based ranking

Unified entry point for TDNN weights. Use `method = "exact"` for the
exact finite-population rank-weight rule
([`weights_tdnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn_exact.md))
and `method = "asymp"` for the asymptotic full-rank approximation
([`weights_tdnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn_asymp.md)).

## Usage

``` r
weights_tdnn(
  X,
  X0 = NULL,
  x0 = NULL,
  distance = "l2",
  s1,
  s2,
  dim = NULL,
  method = c("exact", "asymp"),
  normalize = TRUE,
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- X:

  Numeric matrix (n x p).

- X0:

  Numeric matrix (m x p); mutually exclusive with `x0`.

- x0:

  Numeric length-p single target (convenience for m = 1).

- distance:

  Distance spec or registered name (e.g., `"l2"`).

- s1:

  Positive integer; smaller U-statistic order (`s1 < s2 <= n`).

- s2:

  Positive integer; larger U-statistic order.

- dim:

  Positive integer; localization dimension. Defaults to `ncol(X)` when
  `NULL`.

- method:

  One of `"exact"` or `"asymp"`.

- normalize:

  Logical; renormalize rank table columns to sum 1 (default `TRUE`).

- ctx:

  Optional context passed through to distance and rank-weight
  evaluators.

- return:

  Character vector of extras to include; any of
  `c("W", "ord", "ranks", "dx")`.

## Value

A `weights_result` object.

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

## See also

[`weights_tdnn_exact()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn_exact.md),
[`weights_tdnn_asymp()`](https://JakobJuergens.github.io/RdistNN/reference/weights_tdnn_asymp.md)

## Examples

``` r
set.seed(1)
X  <- matrix(rnorm(30), 15, 2)
x0 <- c(0, 0)
out <- weights_tdnn(X, x0 = x0, s1 = 3L, s2 = 6L, method = "exact")
colSums(out$W)  # 1
#> [1] 1
```
