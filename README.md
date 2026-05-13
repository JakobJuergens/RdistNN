
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RdistNN

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/JakobJuergens/RdistNN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JakobJuergens/RdistNN/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/JakobJuergens/RdistNN/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/JakobJuergens/RdistNN/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

**RdistNN** is a modular toolkit for distance-localized estimation with
DNN rank weights and generalized U-statistics.

- Core distance/rank/weight primitives (`weights_dnn()`,
  `weights_matrix()`, custom registries)
- User-level wrappers for localized estimation: `np_dnn()`,
  `cate_dnn()`, `plqr_dnn()`, and `gmm_dnn()`
- Unified diagnostics for fitted local estimators via `diagnostics()`
- Built-in post-estimation methods for CATE (`summary()`, `confint()`,
  `predict()`)
- Fast delete-1 jackknife options (`jackknife = "fast"` and `"loo"`) for
  standard errors
- Simulation infrastructure for the paper’s CATE study:
  `sim_cate_data()`, `sim_cate_one_rep()`, `sim_cate_run()`,
  `summarise_sim_results()`
- Local IV-GMM estimator and simulation layer: `gmm_model_dr_iv()`,
  `gmm_dnn()`, `sim_iv_data()`, `sim_iv_run()`, `iv_demo_data()`

> Status: active research package. APIs may continue to evolve.

For exact DNN weights, the package works with the full rank permutation
at each target. In particular, exact single-scale DNN with `k = 1` has
support through rank `n - s + 1`; it is not an `s`-nearest-neighbor
truncation. The implementation also supports fixed top-`k` DNN kernels.
Extending the fixed-`k` theory from the 1-NN derivation is essentially
direct, but the current accompanying manuscript derivations are written
for `k = 1`; use larger fixed `k` values as implementation-supported
generalizations until that proof text is added.

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("JakobJuergens/RdistNN")
```

## Quick Start

### 1) DNN Weights and Local Means (`np_dnn`)

``` r
library(RdistNN)

set.seed(1)
demo_np <- np_demo_data(n = 250, p = 2, design = "sinusoidal", seed = 1)

fit_np <- np_dnn(
  data = demo_np$data,
  x_cols = c("X1", "X2"),
  y_col = "Y",
  x0 = c(0, 0),
  distance = "l2",
  method = "exact",
  k = 2L, s = 3L,
  se = TRUE,
  jackknife = "fast",
  return = c("W", "ess", "neighbors")
)

fit_np$mu
fit_np$se
diagnostics(fit_np)
```

### 2) CATE Estimation (`cate_dnn`) and Methods

``` r
library(RdistNN)

demo_cate <- cate_demo_data(n = 300, p = 2, design = "hetero", seed = 1)
df <- demo_cate$data
df$m1_oracle <- df$m0_oracle + df$tau_oracle

fit_cate <- cate_dnn(
  data = df,
  x0 = rep(0, 2),
  x_cols = c("X1", "X2"),
  y_col = "Y",
  t_col = "D",
  m0_col = "m0_oracle",
  m1_col = "m1_oracle",
  e_col = "e_oracle",
  distance = "l2",
  method = "exact",
  k = 2L, s = 3L,
  se = TRUE,
  jackknife = "fast",
  return = c("W", "psi", "nuisance")
)

summary(fit_cate)
confint(fit_cate)
diagnostics(fit_cate)
predict(fit_cate, x0 = c(0.2, -0.1), interval = "confidence")
```

### 3) Local PLQR Orthogonal Moments (`plqr_dnn`)

``` r
library(RdistNN)

demo_plqr <- plqr_demo_data(n = 250, p = 2, p_d = 2, tau = 0.5, seed = 1)

fit_plqr <- plqr_dnn(
  data = demo_plqr$data,
  x_cols = c("X1", "X2"),
  D_cols = c("D1", "D2"),
  q_col = "q_oracle",
  mD_cols = c("mD1", "mD2"),
  y_col = "Y",
  beta = demo_plqr$truth$beta_true,
  tau = 0.5,
  x0 = c(0, 0),
  distance = "l2",
  method = "exact",
  k = 2L, s = 3L,
  se = TRUE,
  jackknife = "fast",
  return = c("W", "g", "ess")
)

fit_plqr$moment
fit_plqr$se
diagnostics(fit_plqr)
```

### 4) CATE Simulation Study

``` r
library(RdistNN)

# Draw one dataset from a named DGP (oracle nuisances included)
dgp <- sim_cate_data("dgp2_nonlinear", n = 300, kx = 2, seed = 1)
dgp$tau_x0   # true CATE at x0 = (0, 0)

# Run one replicate and get bias / SE / coverage metrics
r <- sim_cate_one_rep(n = 300, kx = 2, k = 1, s = 12,
                      design = "dgp2_nonlinear", seed = 1)
r

# Small Monte Carlo grid (use n_reps = 500 for paper results)
raw <- sim_cate_run(
  designs = c("dgp1_linear", "dgp2_nonlinear"),
  n_vals  = 300L,
  n_reps  = 5L,
  verbose = FALSE
)
summarise_sim_results(raw)
```

### 5) Local IV-GMM Estimation and Simulation

``` r
library(RdistNN)

# Draw one local IV-GMM dataset with two valid instruments
iv <- sim_iv_data(
  "dgp2_iv_homog_overid",
  n = 300,
  kx = 2,
  q = 2,
  seed = 1
)

names(iv$data)
iv$theta_x0

model_iv <- gmm_model_dr_iv(
  z_cols = c("Z1", "Z2"),
  pi_cols = c("piZ1_oracle", "piZ2_oracle"),
  muY_col = "muY_oracle",
  muD_col = "muD_oracle"
)

fit_iv <- gmm_dnn(
  data = iv$data,
  x_cols = c("X1", "X2"),
  x0 = iv$x0,
  model = model_iv,
  k = 2L, s = 12L,
  return = "ess"
)

fit_iv$theta
fit_iv$J

# Smaller demo-data wrapper for examples and tests
demo_iv <- iv_demo_data(n = 120, p = 2, q = 2, design = "homog", seed = 1)
head(demo_iv$data)

# Small Monte Carlo grid
raw_iv <- sim_iv_run(
  designs = "dgp2_iv_homog_overid",
  n_vals = 300L,
  q_vals = 2L,
  n_reps = 5L,
  methods = c("dnn_iv_oracle", "iv_oracle"),
  verbose = FALSE
)
summarise_iv_results(raw_iv)
```

### 6) Custom Distances, Rank Weights, and Scores

``` r
library(RdistNN)

register_distance("linf_demo", function(X, x0, ctx = NULL) {
  apply(abs(X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)), 1L, max)
})

register_rank("inv_rank_demo", function(ranks) {
  w <- 1 / as.numeric(ranks)
  w / sum(w)
})

X <- as.matrix(np_demo_data(n = 80, p = 2, seed = 1)$data[, c("X1", "X2")])
W <- weights_matrix(
  X = X,
  x0 = c(0, 0),
  ranking = "linf_demo",
  weight = get_rank("inv_rank_demo"),
  return = "W"
)$W
dim(W)

sc_affine <- new_score(
  name = "affine_order1",
  order = 1L,
  eval_param = function(idxs, data, theta, ctx, rng = NULL) {
    theta[1] * data$Y[idxs] + theta[2]
  }
)

path <- score_diff_path(sc_affine, data = data.frame(Y = c(1, 2, 4)), method = "fd")
path$grad(2L, theta = c(0.5, 1))
```

The full extension guide is
`vignette("extending", package = "RdistNN")`. It covers `setup()` /
`eval_batch()` patterns, custom ranking inputs to `weights_matrix()`,
parameter-aware scores, and debugging with
`return = c("dx", "ord", "ranks")`.

## Documentation Map

- Start with `vignette("getting-started", package = "RdistNN")`
- `vignette("diagnostics-reporting", package = "RdistNN")` for
  fitted-object diagnostics
- `vignette("dnn-weights", package = "RdistNN")` for exact/asymptotic
  DNN weight rules
- `vignette("conditional-z-estimation", package = "RdistNN")` for
  end-to-end estimator workflows
- `vignette("local-iv-estimation", package = "RdistNN")` for local
  IV-GMM estimation and simulation
- `vignette("extending", package = "RdistNN")` for custom
  distance/rank/score components
- `vignette("jackknife-gen-ustats", package = "RdistNN")` for jackknife
  background
- `vignette("paper-cate-simulation", package = "RdistNN")` for the
  paper’s CATE simulation infrastructure

## Development Notes

- Engine selection is automatic: compiled code is used when available,
  otherwise robust R fallbacks are used.
- Deterministic tie handling is built in for reproducibility.
- Tests use `testthat` (edition 3). Run `devtools::test()`.
- Rebuild this README from source with `devtools::build_readme()`.

## Getting Help

Please open an issue with a minimal reproducible example.

------------------------------------------------------------------------

*This package accompanies ongoing work on generalized U-statistics,
jackknife variance estimation, and distance-based localizers.*
