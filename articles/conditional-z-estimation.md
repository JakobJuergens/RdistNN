# Conditional Z-Estimation with DNN weights

## Overview

This article focuses on CATE and nonparametric regression. Start with
[Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md)
for executable examples, or the [CLATE
guide](https://JakobJuergens.github.io/RdistNN/articles/conditional-late.md)
for supplied-prediction binary-IV and paired-signal workflows. This
article gives additional moment details, including **experimental
PLQR**:

- [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md)
  for localized conditional means
- [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
  for conditional average treatment effects
- [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
  for experimental moment evaluation at supplied coefficients

Each wrapper supports:

- matrix and data-frame interfaces
- DNN method choice (`"exact"` or `"asymp"`)
- delete-1 jackknife SEs (`"fast"`, `"loo"`, `"none"`)
- modular kernel-order tuning via `tuning = tuning_dnn_paper(...)`
- fitted-object diagnostics via
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)

The wrappers implement fixed top-`k` DNN kernels for `1 <= k <= s`. The
accompanying manuscript derivations are stated for `k = 1`; larger fixed
`k` is an implementation-supported generalization.

## Local Mean (`np_dnn`)

``` r

demo_np <- np_demo_data(n = 250, p = 2, design = "sinusoidal", seed = 1)

fit_np <- np_dnn(
  data = demo_np$data,
  x_cols = c("X1", "X2"),
  y_col = "Y",
  x0 = c(0, 0),
  method = "exact",
  k = 1L,
  tuning = tuning_default,
  se = TRUE,
  jackknife = "fast",
  return = c("W", "ess", "tuning")
)

fit_np$mu
#> [1] -0.02251305
fit_np$se
#> [1] 0.1782743
diagnostics(fit_np)
#> <dnn_diagnostics>
#>   targets: 1
#>   n_eff  : min=18.6582, median=18.6582, max=18.6582
#>   weights: max(top5 share)=0.4122, max(weight)=0.1
#>   dist   : available for 0/1 target(s)
#>   outcome: available for 1/1 target(s)
```

## CATE (`cate_dnn`)

[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
supports either:

- Full AIPW nuisance inputs (`m0`, `m1`, `e`)
- Robinson-form nuisance inputs (`m`, `e`)
- Direct `pseudo_outcome`

``` r

demo_cate <- cate_demo_data(n = 300, p = 2, design = "hetero", seed = 1)
df <- demo_cate$data
df$m1_oracle <- df$m0_oracle + df$tau_oracle

fit_cate <- cate_dnn(
  data = df,
  x_cols = c("X1", "X2"),
  y_col = "Y",
  t_col = "D",
  m0_col = "m0_oracle",
  m1_col = "m1_oracle",
  e_col = "e_oracle",
  x0 = c(0, 0),
  method = "exact",
  k = 1L,
  tuning = tuning_default,
  se = TRUE,
  return = c("W", "psi", "nuisance", "tuning")
)

summary(fit_cate)
#> <summary.cate_result>
#>   targets: 1
#>          estimate std.error statistic p.value conf.low conf.high
#> target_1   0.5512    0.3181     1.733 0.08313 -0.07227     1.175
confint(fit_cate)
#>             conf.low conf.high
#> target_1 -0.07226519  1.174624
diagnostics(fit_cate)
#> <dnn_diagnostics>
#>   targets: 1
#>   n_eff  : min=20.098, median=20.098, max=20.098
#>   weights: max(top5 share)=0.3894, max(weight)=0.0933
#>   dist   : available for 0/1 target(s)
#>   overlap: available for 1/1 target(s)
predict(fit_cate, x0 = c(0.2, -0.1), interval = "confidence")
#>            fit       lwr      upr
#> [1,] 0.8336633 0.3016249 1.365702
```

## Experimental PLQR moments and coefficient uncertainty (`plqr_dnn`)

All PLQR functions, including score and demo helpers, methods, and
specialized plots, are experimental. Their APIs and inference interfaces
may change. The current wrapper does not estimate a local coefficient
root.

``` r

demo_plqr <- plqr_demo_data(n = 250, p = 2, p_d = 2, tau = 0.5, seed = 1)
# demo_plqr$data$density_oracle holds the true threshold density of the
# demo DGP (constant, 4 * tau * (1 - tau) / sigma).

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
  method = "exact",
  k = 1L,
  tuning = tuning_default,
  se = TRUE,
  return = c("W", "score", "ess", "tuning")
)

summary(fit_plqr)
confint(fit_plqr)
predict(fit_plqr, x0 = c(0.25, -0.25), interval = "confidence")
fit_plqr$moment
diagnostics(fit_plqr)
```

### PLQR nuisances and inference

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
localizes the moment

``` math
\mathbb{E}\!\left[
  \left(\tau - \mathbf{1}\{Y \le D^\top \theta + g_0(U)\}\right)
  \left(D - m_D(U)\right)
  \;\middle|\; X \approx x_0
\right]
```

at the supplied coefficient value `theta_eval`. Here `g` is the
partially-linear model’s nonparametric shape nuisance (paper notation
$`g_0(U)`$), not itself a quantile — the sum $`D'\theta_0 + g_0`$ is the
conditional quantile. `mD` residualizes the design
($`m_0(U) = E[D \mid U]`$). Both are supplied by the user; no nuisance
model is fit internally, and no root solving is performed: the result’s
`theta` repeats the supplied `theta_eval`
(`meta$theta_source = "supplied"`).

Coefficient-level standard errors and confidence intervals require
supplying `density`/`density_col` (an estimate of the local conditional
density of `Y` at the fitted quantile). This is not an optional
convenience: every valid coefficient route needs the density-backed
sandwich bread. Two coefficient-meat routes are available through
`theta_meat`: the analytic score covariance (default) and a full
delete-1 jackknife covariance of the localized vector moment
(`theta_meat = "jackknife"`, stored in `moment_vcov`). A fit without a
density input can still report `moment_se` (a componentwise jackknife SE
on the *moment*, useful as a diagnostic) but
[`confint()`](https://rdrr.io/r/stats/confint.html)/`predict(se.fit = TRUE)`
will error, by design, until refit with a density input — a
componentwise `moment_se` alone is never coefficient inference.
Coefficient-level results live in `theta_se` / `theta_vcov`; the routes
are never mixed in one field.

## Diagnostics and Reporting

The three wrappers return enough retained state for
[`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
to compute fit-level summaries without rerunning the estimator.

``` r

d_np <- diagnostics(fit_np)
d_cate <- diagnostics(fit_cate)
d_plqr <- diagnostics(fit_plqr)

d_np[, c("target_id", "n_eff", "weighted_mean_y")]
d_cate[, c("target_id", "treated_mass", "control_mass", "treated_share")]
d_plqr[, c("target_id", "gram_min_eig", "gram_condition_number")]
```

General diagnostics summarize weight concentration, realized distance
radii when `dx` was retained, and tuning metadata. Estimator-specific
diagnostics add NP outcome summaries, CATE overlap summaries, and PLQR
residual-regressor Gram conditioning.

## Notes on Standard Errors

[`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md)
and
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
always use delete-1 jackknife SEs when `se = TRUE`.

- `"fast"` reuses full-sample ranks and applies rank updates per
  deletion.
- `"loo"` recomputes leave-one-out weights directly
  (reference-equivalent).
- `"none"` skips SE computation.

`s < n` is required when jackknife SEs are requested.

Experimental
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
differs: with `density`/`density_col`, `se = TRUE` reports a
density-backed coefficient sandwich. Its meat uses either the analytic
score covariance (default) or the vector-moment jackknife covariance
with `theta_meat = "jackknife"`. Without density, `moment_se` is a
moment jackknife diagnostic and does not support coefficient intervals.
See “PLQR nuisances and inference” above.

## Related Articles

- [Notation and package
  conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md):
  [`vignette("notation-and-conventions", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md)
  maps arguments, paper symbols, fitted-object fields, and uncertainty
  scales.
- [Diagnostics and
  reporting](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md):
  [`vignette("diagnostics-reporting", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md)
  documents the diagnostics available for each estimator.
- [Kernel-order
  tuning](https://JakobJuergens.github.io/RdistNN/articles/kernel-order-tuning.md):
  [`vignette("kernel-order-tuning", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/kernel-order-tuning.md)
  explains how `s` is selected when it is not supplied directly.
