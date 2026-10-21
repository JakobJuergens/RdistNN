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
- for the exact finite-population rule, support in rank order is
  `j <= n - s + k`, and in particular `j <= n - s + 1` when `k = 1`;
- the asymptotic Beta-increment approximation is typically treated as a
  full-support rank-weight sequence.

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
