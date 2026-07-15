# DNN weights: exact and asymptotic L-statistics

## Overview

`RdistNN` provides two DNN rank-weight constructions:

- `method = "exact"`: finite-sample exact DNN coefficients
- `method = "asymp"`: asymptotic DNN coefficients

Use
[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md)
to switch between them with a common interface.

Two support facts matter in practice:

- exact DNN weights are computed from the **full rank permutation** for
  each target, not by keeping only the first `s` observations;
- for the exact finite-population rule, support in rank order is
  `j <= n - s + k`, and in particular `j <= n - s + 1` when `k = 1`;
- the asymptotic Beta-increment approximation is typically treated as a
  full-support rank-weight sequence.

Fixed top-`k` DNN kernels are fully implemented for `1 <= k <= s`. The
fixed-`k` extension of the 1-NN theory is essentially direct, but the
accompanying manuscript currently writes out the theoretical derivations
only for `k = 1`.

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
colSums(wt_asymp$W)
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
```

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
```

For moderate samples and small `s`, exact weights are often preferred
when you want finite-sample fidelity. Asymptotic weights can be cheaper
and are often a practical default for larger-scale usage.
