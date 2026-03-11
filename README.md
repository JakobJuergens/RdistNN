
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RdistNN

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
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
n <- 500
x <- runif(n)
y <- sin(2 * pi * x) + rnorm(n, sd = 0.2)
X <- matrix(x, ncol = 1)
x0 <- 0.4

# Exact DNN weights at x0 using k = 1, s = 1
w1 <- weights_dnn_exact(X, x0 = x0, distance = "l2", k = 1L, s = 1L)$W[, 1]

# Or use k = 5 with subset order s = 5
wk <- weights_dnn_exact(X, x0 = x0, distance = "l2", k = 5L, s = 5L)$W[, 1]

# Local average at x0
mu_hat <- sum(wk * y)
mu_hat
```

### 2) PLQR score components (work in progress)

``` r
library(RdistNN)

demo <- plqr_demo_data(n = 200, p = 2, p_d = 2, tau = 0.5, seed = 1)
sc <- score_plqr_component(
  x0 = rep(0, 2),
  x_cols = c("X1", "X2"),
  D_cols = c("D1", "D2"),
  q_col = "q_oracle",
  mD_cols = c("mD1", "mD2"),
  beta = demo$truth$beta_true,
  j = 1L,
  tau = 0.5,
  s = 3L,
  k = 1L
)
ctx <- score_setup(sc, demo$data)
subs <- plqr_demo_subsets(nrow(demo$data), K = 5, s = 3L, seed = 1)
score_eval_batch(sc, subs, demo$data, ctx)
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
