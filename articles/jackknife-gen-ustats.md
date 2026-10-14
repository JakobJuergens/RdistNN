# Jackknife for Hájek-dominated generalized U-statistics (overview)

## Overview

`RdistNN` includes two delete-1 jackknife utilities for rank-weighted
functionals:

- [`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md):
  optimized implementation (compiled when available)
- [`jackknife_slow_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_slow_rank.md):
  reference implementation for validation

Both return:

- `theta0`: baseline aggregate
- `theta_loo`: all leave-one-out aggregates
- `info`: metadata (`s`, number of targets, output dimension,
  normalization)

Important distinction:

- [`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md)
  and
  [`jackknife_slow_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_slow_rank.md)
  are **generic** utilities for compact-support rank-weighted L-stat
  rules;
- the user-level exact DNN wrappers do **not** truncate to the first `s`
  observations. They keep the full rank matrix and reevaluate the DNN
  rank-weight rule after delete-one updates.

## Example Setup

``` r

n <- 40
p <- 2
s <- 3L

X <- matrix(rnorm(n * p), n, p)
Y <- rnorm(n)

# Toy compact-support example: build only the first s+1 ranked observations
# for each target point. This illustrates the generic L-stat utilities below;
# it is not the internal exact-DNN jackknife path used by np_dnn()/cate_dnn()/plqr_dnn().
wt <- weights_dnn(
  X = X, X0 = X,
  method = "exact",
  distance = "l2",
  k = s, s = s,
  return = "ord"
)
ranks <- t(wt$ord[seq_len(s + 1L), , drop = FALSE])  # T x (s+1)

# L-stat coefficients
a <- rep(1, s)

# Scalar signal represented as n x 1 matrix
phi <- cbind(Y)
```

## Fast vs Slow Jackknife

``` r

jk_fast <- jackknife_fast_rank(ranks = ranks, a = a, phi = phi, normalize = TRUE)
jk_slow <- jackknife_slow_rank(ranks = ranks, a = a, phi = phi, normalize = TRUE)

all.equal(jk_fast$theta0, jk_slow$theta0, tolerance = 1e-12)
all.equal(jk_fast$theta_loo, jk_slow$theta_loo, tolerance = 1e-12)
```

## Turning LOO Aggregates into Standard Errors

``` r

theta_loo <- jk_fast$theta_loo[, 1]
theta_bar <- mean(theta_loo)
n_ref <- length(theta_loo)

var_jk <- ((n_ref - 1) / n_ref) * sum((theta_loo - theta_bar)^2)
se_jk <- sqrt(var_jk)
se_jk
```

## Notes

- [`jackknife_fast_rank()`](https://JakobJuergens.github.io/RdistNN/reference/jackknife_fast_rank.md)
  supports both matrix `phi` and function-valued `phi(j, t)`.
- The slow version is useful for debugging and parity checks in tests.
- For user-level wrappers
  ([`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)),
  these steps are handled internally when `se = TRUE`.

## Related Articles

- [Notation and package
  conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md):
  [`vignette("notation-and-conventions", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md)
  distinguishes DNN kernel order, localized moments, and the uncertainty
  fields returned by user-level estimators.
- [DNN
  weights](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md):
  [`vignette("dnn-weights", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md)
  describes the exact rank weights used by the jackknife path.
- [Conditional
  Z-estimation](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md):
  [`vignette("conditional-z-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md)
  shows the high-level wrappers that invoke these routines internally.
