# RdistNN

**RdistNN** implements distributional nearest-neighbor (DNN)
localization for local conditional moment and Z-estimation problems. The
pre-public interface overhaul is implemented; the package now uses
consistent estimator, tuning, diagnostic, plotting, and simulation
contracts across its primary workflows.

- Core distance/rank/weight primitives
  ([`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
  [`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md),
  custom registries)
- User-level wrappers for localized estimation:
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  and
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
- Unified diagnostics for fitted local estimators via
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
- Built-in post-estimation methods
  ([`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html))
- Tuning-grid inspection through
  [`candidate_fit()`](https://JakobJuergens.github.io/RdistNN/reference/candidate_fit.md),
  [`candidate_fits()`](https://JakobJuergens.github.io/RdistNN/reference/candidate_fit.md),
  and
  [`plot_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/plot_tuning.md)
- ggplot2 postestimation views for estimates, overlap, weights,
  influence, conditioning, calibration, and dashboards
- Fast delete-1 jackknife options (`jackknife = "fast"` and `"loo"`) for
  standard errors
- Simulation infrastructure for the paper’s CATE study:
  [`sim_cate_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_data.md),
  [`sim_cate_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_one_rep.md),
  [`sim_cate_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_run.md),
  [`summarise_sim_results()`](https://JakobJuergens.github.io/RdistNN/reference/summarise_sim_results.md)
- Experimental companion functionality: local IV-GMM wrappers and
  simulation helpers
  ([`gmm_model_dr_iv()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_model_dr_iv.md),
  [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md),
  [`sim_iv_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_data.md),
  [`sim_iv_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_run.md),
  [`iv_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/iv_demo_data.md))
  are included for ongoing companion-paper development; their API and
  theorem alignment may change before the GMM paper is released

> Status: active research package. APIs may continue to evolve.

This package accompanies [Orthogonal Inference for Conditional
Z-Estimation without
Cross-Fitting](https://github.com/JakobJuergens/Cond_Z_Estimation).

## API tiers

Not all exported functions carry the same stability promise:

| Tier | Meaning | Examples | Public docs |
|----|----|----|----|
| Stable user API | expected to survive the public release | [`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md), [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md), [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md), [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md), fitted-object methods, diagnostics, tuning, and plotting helpers | README + main reference groups |
| Developer API | extension points for custom components | distance / rank-weight / score specs and registries | “Developer API” reference groups |
| Paper replication API | tied to the paper’s simulation study | `sim_cate_*`, [`summarise_sim_results()`](https://JakobJuergens.github.io/RdistNN/reference/summarise_sim_results.md) | replication reference group + vignette |
| Experimental companion API | may change before the companion papers are released | [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md), `sim_iv_*`, TDNN/jackknife infrastructure | “Experimental”/“Advanced” reference groups |
| Internal | not exported | low-level distance/order/rank-table plumbing | none |

For exact DNN weights, the package works with the full rank permutation
at each target. Exact single-scale DNN has finite-population support
through rank `n - s + k`; in the paper-default `k = 1` case this is
`n - s + 1`. Thus `s` is a subsample/localization order, not an
`s`-nearest-neighbor truncation. The optional `neighbors` extra reports
the first `k` nearest observations used inside each size-`s` subsample,
not the full DNN weight support. The implementation also supports fixed
top-`k` DNN kernels. Extending the fixed-`k` theory from the 1-NN
derivation is essentially direct, but the current accompanying
manuscript derivations are written for `k = 1`; use larger fixed `k`
values as implementation-supported generalizations until that proof text
is added.

## Installation

Install the development version from GitHub:

``` r

# install.packages("pak")
pak::pak("JakobJuergens/RdistNN")
```

## Platform compatibility

R package binaries install without a local compiler when a matching
binary is available. Installation from source requires a working C++17
toolchain and `Rcpp` on Windows, macOS, and Linux. OpenMP is optional:
the package uses R’s portable `SHLIB_OPENMP_CXXFLAGS`, so compiled
kernels run serially on toolchains where R does not provide OpenMP
flags.

The package includes R implementations for supported distance, score,
and jackknife operations and uses them when the corresponding compiled
fast path is unavailable at runtime. Set
`options(RdistNN.use_compiled = FALSE)` to force those R paths for
reproducibility or diagnosis; CI checks compiled/fallback parity on
every supported platform. Those algorithmic fallbacks do not turn the
package into a no-compiler source package: the package itself contains
registered C++ code and must still be compiled when installed from
source.

The full GitHub Actions release matrix covers Linux x86-64, Windows
x86-64, macOS arm64, and macOS x86-64, plus R release, devel, and oldrel
on Linux. Pull requests, pushes to the default branch, and manual runs
all execute the full matrix so platform regressions are blocked before
merge.

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
  # Toy choice for a tiny example; use tuning_dnn_paper() in paper-style runs.
  k = 1L, s = 3L,
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
  # Toy choice for a tiny example; use tuning_dnn_paper() in paper-style runs.
  k = 1L, s = 3L,
  se = TRUE,
  jackknife = "fast",
  return = c("W", "psi", "nuisance")
)

summary(fit_cate)
confint(fit_cate)
diagnostics(fit_cate)
predict(fit_cate, x0 = c(0.2, -0.1), interval = "confidence")
```

#### CATE inference contract

Neyman orthogonality controls first-order population nuisance drift. It
does not, by itself, validate arbitrary nuisance predictions estimated
on the same observations used by
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md).
Valid full-sample reuse additionally requires the nuisance-rate and
full-state stability conditions in the paper.
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
computes the estimator and requested jackknife procedure but cannot
diagnose whether an arbitrary learner satisfies those conditions.
Supported inputs include oracle quantities, cross-fitted or out-of-fold
predictions, and full-sample predictions from learners for which the
paper’s rate and stability conditions are credible.

### 3) Local PLQR Orthogonal Moments (`plqr_dnn`)

PLQR support currently evaluates localized orthogonal PLQR moments and
direct-sandwich coefficient uncertainty at supplied coefficient values.
Local root-solving for target-specific PLQR coefficients is under
development. Coefficient-level standard errors and confidence intervals
require threshold-density values via `density` / `density_col`.

``` r

library(RdistNN)

demo_plqr <- plqr_demo_data(n = 250, p = 2, p_d = 2, tau = 0.5, seed = 1)
# demo_plqr$data already includes density_oracle, the true (constant)
# conditional density of Y at the quantile threshold for the demo DGP.

fit_plqr <- plqr_dnn(
  data = demo_plqr$data,
  x_cols = c("X1", "X2"),
  D_cols = c("D1", "D2"),
  g_col = "g_oracle",
  mD_cols = c("mD1", "mD2"),
  y_col = "Y",
  density_col = "density_oracle",
  theta_eval = demo_plqr$truth$theta_true,
  tau = 0.5,
  x0 = c(0, 0),
  distance = "l2",
  method = "exact",
  k = 1L,
  tuning = tuning_dnn_paper(),
  se = TRUE,
  return = c("W", "score", "ess")
)

summary(fit_plqr)
confint(fit_plqr)
predict(fit_plqr, x0 = c(0.25, -0.25), interval = "confidence")
diagnostics(fit_plqr)
```

[`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md)
uses the paper’s finite-sample default `lambda = 0.4` and sensitivity
grid `c(0.3, 0.4, 0.5, 0.75)`, with
`s(lambda) = floor(lambda * n^((kx + 1)/(kx + 2)))`. Explicit `lambda`,
`lambda_grid`, fixed `s`, and custom tuning specifications remain
supported; `lambda = 1` is a more aggressive localization benchmark.

Without a density input,
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
still evaluates the localized moment and reports `moment_se`, a
jackknife **moment diagnostic**; that is not coefficient-level
inference.

### 4) CATE Simulation Study

``` r

library(RdistNN)

# Draw one dataset from a named DGP (oracle AIPW and nuisance columns included)
dgp <- sim_cate_data("dgp2_nonlinear", n = 300, kx = 2, seed = 1)
dgp$tau_x0   # true CATE at the default x0 = (1, 0)

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

For DNN and local-constant rows, the raw and summary outputs also expose
the realized DNN-weighted truth, realized localization drift, estimator
remainder, physical neighborhood radius, arm-specific effective sample
sizes, and studentized finite-sample shape. The weighted truth
`sum_i W_i,s(x0) tau_0(X_i)` is random and is not the deterministic
population branch `alpha_s(x0)`. With `x0` and `(n, s)` fixed,
expectation over the i.i.d. covariate sample and complete-DNN
exchangeability give `E[sum_i W_i,s(x0) tau_0(X_i)] = alpha_s(x0)`. The
estimator remainder combines centered score variation with fitted-state
effects; paired feasible/oracle comparisons isolate the latter. These
fields keep localization, score variation, and nuisance fitting distinct
instead of treating coverage as a single undifferentiated diagnostic.

### 5) Local IV-GMM (experimental)

Experimental local IV-GMM wrappers
([`gmm_model_dr_iv()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_model_dr_iv.md),
[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md))
and simulation helpers (`sim_iv_*`,
[`iv_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/iv_demo_data.md))
are included for companion-paper development. Their API may change
before the GMM paper is released. See
[`vignette("local-iv-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/local-iv-estimation.md)
for the full workflow.

### 6) Custom Distances, Rank Weights, and Scores

``` r

library(RdistNN)

register_distance("linf_demo", function(X, x0, ctx = NULL) {
  apply(abs(X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)), 1L, max)
})

register_rank_weight("inv_rank_demo", function(ranks) {
  w <- 1 / as.numeric(ranks)
  w / sum(w)
})

X <- as.matrix(np_demo_data(n = 80, p = 2, seed = 1)$data[, c("X1", "X2")])
W <- weights_matrix(
  X = X,
  x0 = c(0, 0),
  ranking = "linf_demo",
  weight = get_rank_weight("inv_rank_demo"),
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
[`vignette("extending", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/extending.md).
It covers `setup()` / `eval_batch()` patterns, custom ranking inputs to
[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md),
parameter-aware scores, and debugging with
`return = c("dx", "ord", "ranks")`.

## Visualisation

Fitted objects have ggplot2-based
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)/`autoplot()`
methods (ggplot2 is a suggested dependency). The flagship view is the
CATE profile: point estimates with jackknife confidence bounds along a
covariate profile, a zero reference line, and an honesty caption
reporting the effective local sample size and the kernel-order
provenance.

``` r

library(RdistNN)

set.seed(42)
n <- 1000
X <- cbind(X1 = runif(n, -2, 2), X2 = rnorm(n))
e <- plogis(0.4 * X[, "X1"])                    # oracle propensity
tau <- 1 + sin(X[, "X1"])                       # true CATE along X1
m0 <- X[, "X1"] - 0.5 * X[, "X2"]               # oracle control mean
treat <- as.integer(runif(n) < e)
Y <- m0 + treat * tau + rnorm(n, sd = 0.5)

df <- data.frame(X, T = treat, Y = Y, e = e, m0 = m0, m1 = m0 + tau)
X0 <- cbind(X1 = seq(-1.8, 1.8, length.out = 21), X2 = 0)

fit <- cate_dnn(
  data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
  e_col = "e", m0_col = "m0", m1_col = "m1",
  X0 = X0, k = 1L, se = TRUE, jackknife = "fast"
)

plot(fit)
```

![](reference/figures/README-cate-profile-1.png)

The plotting surface mirrors the fitted-object diagnostic contract:
[`plot_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_estimates.md)
handles coefficient profiles and intervals,
[`plot_overlap()`](https://JakobJuergens.github.io/RdistNN/reference/plot_overlap.md)
shows local treatment and propensity support,
[`plot_influence()`](https://JakobJuergens.github.io/RdistNN/reference/plot_influence.md)
identifies high-contribution observations,
[`plot_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/plot_tuning.md)
compares feasible tuning candidates, and
[`plot_dashboard()`](https://JakobJuergens.github.io/RdistNN/reference/plot_dashboard.md)
assembles the main review panels.
[`plot_conditioning()`](https://JakobJuergens.github.io/RdistNN/reference/plot_conditioning.md)
and
[`plot_calibration()`](https://JakobJuergens.github.io/RdistNN/reference/plot_calibration.md)
provide specialized PLQR/GMM views.

Weight-geometry views show what the estimator averaged over:
`plot_weights(fit, type = "kernel")` draws the effective data-driven
kernel (weight against distance), `type = "profile"` the decay by rank,
and `type = "cumulative"` the concentration of the weight mass. All
panels are regular ggplot objects and can be re-themed or extended with
`+` as usual;
[`plot_data_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_data_estimates.md),
[`plot_data_weights()`](https://JakobJuergens.github.io/RdistNN/reference/plot_data_weights.md),
and
[`plot_data_diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/plot_data_diagnostics.md)
expose the underlying tidy frames for custom figures.

## Documentation Map

- Start with
  [`vignette("getting-started", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md)
- [`vignette("diagnostics-reporting", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md)
  for fitted-object diagnostics
- [`vignette("dnn-weights", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md)
  for exact/asymptotic DNN weight rules
- [`vignette("conditional-z-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md)
  for end-to-end estimator workflows
- [`vignette("local-iv-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/local-iv-estimation.md)
  for local IV-GMM estimation and simulation
- [`vignette("extending", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/extending.md)
  for custom distance/rank/score components
- [`vignette("jackknife-gen-ustats", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/jackknife-gen-ustats.md)
  for jackknife background
- [`vignette("paper-cate-simulation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/paper-cate-simulation.md)
  for the paper’s CATE simulation infrastructure

## Development Notes

- Engine selection is automatic: compiled code is used when available,
  otherwise robust R fallbacks are used.
- Deterministic tie handling is built in for reproducibility.
- Tests use `testthat` (edition 3). Run `devtools::test()` for the suite
  and `devtools::check()` for the release gate.
- The remaining pre-public documentation item is a dedicated notation
  and package-conventions article; the general conditional-Z front end
  and package-level SIPP workflow are longer-term work.
- Rebuild this README from source with `devtools::build_readme()`.

## Getting Help

Please open an issue with a minimal reproducible example.

------------------------------------------------------------------------

*This package accompanies ongoing work on generalized U-statistics,
jackknife variance estimation, and distance-based localizers.*
