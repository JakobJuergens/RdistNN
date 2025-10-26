
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RdistNN

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/JakobJuergens/RdistNN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JakobJuergens/RdistNN/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/JakobJuergens/RdistNN/graph/badge.svg)](https://app.codecov.io/gh/JakobJuergens/RdistNN)
<!-- badges: end -->

**RdistNN** is a small, focused toolkit for distance-based
nearest-neighbor (DNN) kernels and U-statistic building blocks used in
conditional Z-estimation and partially linear quantile regression
(PLQR).

- Fast C++ implementations (with automatic R fallbacks) for core DNN
  primitives  
- Deterministic tie handling and permutation-invariance checks  
- Score factories for PLQR and CATE-style estimators, with tests and CI

> Status: early research package. APIs may change.

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("JakobJuergens/RdistNN")
```

## Quick start

Non-evaluating examples that show intended usage without requiring all
features to be present yet.

### 1) Local mean with DNN weights

``` r
library(RdistNN)

set.seed(1)
n  <- 500
x  <- runif(n)
y  <- sin(2*pi*x) + rnorm(n, sd = 0.2)
X  <- matrix(x, ncol = 1)
x0 <- 0.4

# Nearest-neighbor weight vector at x0 using k = 1 (DNN-1)
w1 <- dnn1(x0, X)           # length n, nonnegative, sums to 1

# Or k-NN weights (k ≥ 1); when k = 1 this should equal dnn1()
wk <- dnnk(x0, X, k = 5)

# Local average at x0
mu_hat <- sum(wk * y)
mu_hat
```

### 2) PLQR score components (work in progress)

``` r
library(RdistNN)

# Suppose you have (Y, D, X, Z), a target quantile tau, and a component index j.
tau <- 0.5
j   <- 1L

# Factory returns a function you can map over subsets / batches
sc <- score_plqr_component(tau = tau, j = j)

# Minimal illustration:
# sc(data = list(Y = Y, D = D, X = X, Z = Z), beta = beta_init)
```

## Design notes

- **Engine selection.** Uses compiled code when available; otherwise
  falls back to R. Unit tests check parity.  

- **Determinism.** Tie handling and ordering are fixed so runs are
  reproducible.  

- **Testing.** `testthat` (edition 3) and GitHub Actions. Run:

  ``` r
  devtools::test()
  ```

- **Docs.** Rebuild this file with:

  ``` r
  devtools::build_readme()
  ```

## Roadmap (short)

- DNN / DNN-k kernels with auto-bandwidth helpers  
- Score factories for CATE and PLQR  
- Jackknife / IJ helpers for generalized U-statistics  
- Vignettes with small simulations once the API stabilizes

## Getting help

Please open an issue with a minimal reproducible example once the public
API stabilizes. Until then, expect frequent changes.

------------------------------------------------------------------------

*This package accompanies ongoing work on generalized U-statistics,
jackknife variance estimation, and distance-based localizers.*
