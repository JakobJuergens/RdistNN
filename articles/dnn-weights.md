# DNN weights and ranks

## Overview

`RdistNN` provides two DNN rank-weight constructions:

- `method = "exact"`: finite-sample exact DNN coefficients
- `method = "asymp"`: asymptotic DNN coefficients

Use
[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md)
to switch between them with a common interface.

The following support conventions matter in practice:

- exact DNN weights are computed from the **full rank permutation** for
  each target, not by keeping only the first `s` observations;
- for the unaveraged exact finite-population rank table, support is
  `j <= n - s + k`, and in particular `j <= n - s + 1` when `k = 1`;
- the asymptotic Beta-increment approximation is typically treated as a
  full-support rank-weight sequence.

At each target, observations with exactly equal represented distance
under the selected metric share the average of their assigned rank
weights. A group crossing the exact table’s support edge can therefore
give observations beyond that rank positive weight. The group’s signed
or unnormalized mass is preserved; individual weights and ESS can
change. No tolerance or rounding groups nearby distances. Distinct
covariate profiles can have the same distance.

Fixed top-`k` DNN kernels are implemented for `1 <= k <= s`. The
accompanying
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
gives the theory for `k = 1`; accepting another `k` does not establish
its inferential justification.

## Basic Usage

``` r

n <- 150
p <- 2
X <- matrix(rnorm(n * p), n, p)
X0 <- matrix(rnorm(8), 4, p)

wt_exact <- weights_dnn(
  X = X, X0 = X0,
  distance = "l2",
  method = "exact",
  k = 1L, s = 4L,
  return = c("W", "ranks", "ord", "dx")
)

wt_asymp <- weights_dnn(
  X = X, X0 = X0,
  distance = "l2",
  method = "asymp",
  k = 1L, s = 4L,
  normalize = TRUE,
  return = "W"
)

colSums(wt_exact$W)
#> [1] 1 1 1 1
colSums(wt_asymp$W)
#> [1] 1 1 1 1
```

## Single-target Convenience

Use `x0` for one target instead of a one-row `X0`.

``` r

x0 <- c(0, 0)
wt_single <- weights_dnn(
  X = X, x0 = x0,
  method = "exact",
  k = 1L, s = 4L,
  return = c("W", "ranks")
)

dim(wt_single$W)  # n x 1
#> [1] 150   1
```

## Exact distance ties

Grouping uses the actual represented output of the selected metric
separately for each target. It averages assigned weights over
observations, including zero entries beyond nominal support; it does not
apply the rank-weight function at an average rank or weight each
distinct covariate profile equally. Returned ranks sort by distance then
original row index. Changing row order only changes this bookkeeping
within a tie; aligned weights, estimates and covariance are invariant up
to numerical precision.

For `X = (-1, 1, 2, 3, 4)`, target zero, `n = 5`, `s = 3` and `k = 1`,
the exact rank table is `(0.6, 0.3, 0.1, 0, 0)`. The first two
observations have the same distance and receive `(0.6 + 0.3)/2 = 0.45`
each:

``` r

X_tied <- matrix(c(-1, 1, 2, 3, 4), ncol = 1)
y_tied <- c(0, 10, 0, 0, 0)
wt_tied <- weights_dnn(X = X_tied, x0 = 0, k = 1L, s = 3L)
#> Warning: weights_dnn(): Exact distance ties detected; weights are averaged
#> within distance blocks, including jackknife deletions. The accompanying theory
#> assumes a continuous distance distribution and does not establish inference
#> under persistent ties.
drop(wt_tied$W)  # 0.45, 0.45, 0.1, 0, 0
#> [1] 0.45 0.45 0.10 0.00 0.00
sum(drop(wt_tied$W) * y_tied)  # 4.5
#> [1] 4.5
```

For a group at ranks `a:b`, every member receives
`sum(w_n[a:b])/(b-a+1)`. Each jackknife deletion removes the actual
observation before applying the selected `n-1` table to the survivors.
Deleting a member leaves `b-a` survivors at ranks `a:(b-1)`, each
receiving `sum(w_(n-1)[a:(b-1)])/(b-a)`; a deleted singleton contributes
no signal. Both the table and the denominator change. In the example the
five deleted estimates are `(7.5, 0, 5, 5, 5)`, whose mean is `4.5`:

``` r

deleted <- vapply(seq_along(y_tied), function(i) {
  keep <- seq_along(y_tied) != i
  w <- weights_dnn(X = X_tied[keep, , drop = FALSE],
                   x0 = 0, k = 1L, s = 3L)$W
  sum(drop(w) * y_tied[keep])
}, numeric(1))
deleted
#> [1] 7.5 0.0 5.0 5.0 5.0
mean(deleted)
#> [1] 4.5
```

For exact DNN, this is the average over uniform orderings within each
distance tie and uniform subsets of size `s`. Within a given subset, a
boundary tie with `b` included observations and `r` remaining neighbor
slots gives each tied member inclusion probability `r/b`, with the usual
`1/k` averaging. For `k = 2`, one closer observation and two tied
boundary observations receive `(1/2, 1/4, 1/4)` when all three are in
the subset. Including every boundary observation with equal weight would
give a different rule. The same block operation applies to asymptotic
tables and signed TDNN weights, preserving their assigned masses.

If all represented distances to a target are equal and the weights are
normalized, each observation receives `1/n`. With frozen geometry, each
survivor receives `1/(n-1)` after a deletion. Increasing `s` cannot
distinguish observations within such a group. Groups can differ at other
targets and under another metric.

The paper assumes a continuous distance distribution, under which the
earlier and current rules agree almost surely. Calls warn when exact
distance ties occur because inference with persistent ties is not
established. Averaging does not establish a uniform mean-squared-error
improvement.

A rank table alone cannot reconstruct these weights without the
per-target groups. Bare supplied `ranks` or `ord` retain strict-order
semantics; a list that also retains `dx` uses those distances for
grouping while preserving the supplied permutation. Fits cache and
validate the selected policy and groups for prediction, candidate views,
covariance, diagnostics and plots. Refit objects created under earlier
tie policies, including version 0.1.2.

## Working with Returned Diagnostics

`return` can include:

- `"W"`: weight matrix (`n x m`)
- `"ranks"`: rank matrix (`n x m`)
- `"ord"`: order matrix (`n x m`)
- `"dx"`: distance matrix (`n x m`)

These are useful when you want to:

- verify deterministic ordering in tie settings
- debug custom distance behavior
- reuse ranks in downstream routines

## Exact vs Asymptotic

``` r

mu_signal <- rnorm(n)
mu_exact <- as.numeric(crossprod(mu_signal, wt_exact$W))
mu_asymp <- as.numeric(crossprod(mu_signal, wt_asymp$W))

cbind(mu_exact, mu_asymp)
#>        mu_exact   mu_asymp
#> [1,] 0.10286550 0.10246327
#> [2,] 0.11298128 0.11244305
#> [3,] 0.21969143 0.21866562
#> [4,] 0.04561642 0.04615826
```

Exact weights are the package default. Setting `method = "asymp"`
changes the weight construction; it is separate from choosing a compiled
or R calculation engine.

## Numerical limits

Asymptotic tables use host Beta-distribution routines. At the smallest
representable scales, rounding the interval endpoints can erase a tail
mass or give R and Stata different zero versus nonzero support.
Averaging several order statistics and combining signed TDNN tables can
also lose tiny masses. This remains an unresolved numerical limitation.
Extreme outcomes can make such differences visible in estimates and
jackknife covariance.

Renormalizing the table or selecting a compiled deletion engine cannot
recover probability already lost during table construction. Exact
weights use a different finite-population rule; switching to them
changes the method rather than repairing an asymptotic table. The
maintainer reference on [weighted
accumulation](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/19_weighted_accumulation.md)
describes the arithmetic and its limits.

## Related Articles

- [Notation and
  conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md):
  orders, support and effective size.
- [Kernel-order
  tuning](https://JakobJuergens.github.io/RdistNN/articles/kernel-order-tuning.md):
  sensitivity comparisons for `s`.
- [Extending
  RdistNN](https://JakobJuergens.github.io/RdistNN/articles/extending.md):
  custom distances and rank-weight maps.
