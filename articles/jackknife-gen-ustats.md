# Using the rank-jackknife utilities

This advanced guide uses
[`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md)
and its reference
[`jackknife_slow_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_slow_rank.md)
for a compact-support rank-weighted rule. Applied users obtain local
standard errors directly from their fitted estimators; see [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md).
Statistical conditions and theory belong to the
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
and its references.

## Supply neighbor indices and a signal

For each target, supply the nearest `support + 1` observation indices.
The extra neighbor replaces a deleted member. Here the rule averages the
first three neighbors at each of two targets:

``` r

set.seed(4)
n <- 40L
support <- 3L
X <- matrix(rnorm(n * 2), n, 2)
X0 <- rbind(c(0, 0), c(1, 0))
Y <- rnorm(n)

dx <- distance_eval_batch(distance_l2(), X, X0)
ord <- apply(dx, 2L, order)
neighbors <- t(ord[seq_len(support + 1L), , drop = FALSE])
a <- rep(1, support)
phi <- cbind(Y)
```

The argument named `ranks` below contains **observation indices**, with
one row per target. This differs from the `n x m` rank permutations used
by
[`jackknife_rank_full()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_rank_full.md).
The compact-support rule is not the exact DNN weight rule: the
unaveraged exact DNN rank table extends to `n-s+k`, rather than the
first `s` observations. Covariate-aware block averaging can extend
nonzero observation weights beyond that rank.

## Compare fast and reference calculations

``` r

fast <- jackknife_fast_rank(ranks = neighbors, a = a, phi = phi)
slow <- jackknife_slow_rank(ranks = neighbors, a = a, phi = phi)
all.equal(fast$theta_loo, slow$theta_loo, tolerance = 1e-12)
#> [1] TRUE
```

`theta0` averages the rule across the supplied target rows; it is not a
vector of separate local fits. These rank-only utilities do not identify
duplicate covariate groups. High-level covariate-aware estimators
additionally average duplicate-block weights in the full sample and
every deletion, with a warning that this practical modification is
outside the paper’s theory. Deleting a member reduces its block size as
well as the reference-sample size; the survivors use the `n-1` table
divided by the remaining block size. The shared estimator path reuses
prefix/postfix deletion sums with a block correction, rather than
rebuilding a weight matrix for each removed row. See the [worked
deletion
example](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.html#identical-covariate-rows).

Block deletions accumulate the original signals with their actual
survivor weights and retain both parts of each compensated prefix/suffix
sum until the final merge. Averaging a block’s signals first can erase a
small observation that survives deletion of a dominant member. Shared
origin selection also accounts for tiny coefficients and nonzero tiny
table masses: exceptional calculations retain the separately rounded raw
products before division by mass. These safeguards use finite double
precision and cannot recover mass already lost in [asymptotic table
construction](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.html#numerical-limits).

`theta_loo` has one row per deleted reference observation and one column
per signal component.

``` r

deleted <- fast$theta_loo[, 1L]
n_ref <- length(deleted)
c(aggregate = fast$theta0,
  standard_error = sqrt((n_ref - 1) / n_ref * sum((deleted - mean(deleted))^2)))
#>      aggregate standard_error 
#>     0.07850962     0.73102688
```

Use a numeric `n x q` signal matrix to retain the full reference-sample
size. Function-valued `phi(j, t)` supports target-specific
contributions, but the reference size is then inferred as `max(ranks)`.
See
[`?jackknife_fast_rank`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md)
for the full contract and
[`?jackknife_rank_full`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_rank_full.md)
for full-rank weight tables. The [notation
guide](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md)
maps the array conventions.
