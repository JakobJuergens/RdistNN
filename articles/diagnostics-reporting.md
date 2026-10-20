# Diagnostics and Reporting

## Overview

[`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
gives a compact table for fitted DNN-localized estimators. The method is
intentionally fit-object based:

- it uses retained weights rather than refitting;
- it uses retained distances only when `dx` was requested;
- it uses hidden estimator state for NP, CATE, and PLQR diagnostics;
- it marks unavailable diagnostics explicitly instead of reconstructing
  missing inputs from assumptions.

The general fields are available for `np_result`, `cate_result`, and
`plqr_result` objects when weights are retained:

- `n_eff`, `max_weight`, and top-weight shares;
- realized rank radii and weighted distance summaries when `dx` is
  retained;
- tuning metadata when a tuning spec was used.

Estimator-specific fields are added by the corresponding method:

- NP: weighted outcome mean, standard deviation, and outcome quantiles;
  signed TDNN weights do not define a probability distribution, so
  weighted SD and quantiles are `NA` when weights are negative. Weight
  concentration summaries then describe signed weights rather than
  probability mass.
- CATE: treated/control mass and propensity summaries;
- PLQR: residual-regressor Gram eigenvalues and condition number.

Diagnostics require `W`. The default `return = "W"` keeps it. If a fit
was created without weights, refit with `return` including `"W"`.

## Nonparametric Regression

``` r

demo_np <- np_demo_data(n = 120, p = 2, design = "sinusoidal", seed = 1)

fit_np <- np_dnn(
  data = demo_np$data,
  x_cols = c("X1", "X2"),
  y_col = "Y",
  x0 = demo_np$x0,
  method = "exact",
  k = 1L,
  s = 4L,
  se = FALSE,
  return = c("W", "dx", "ess")
)

d_np <- diagnostics(fit_np)
d_np[, c(
  "target_id",
  "n_eff",
  "weighted_mean_y",
  "weighted_sd_y",
  "weighted_y_q50"
)]
#> <dnn_diagnostics>
#>   targets: 1
#>   n_eff  : min=51.8421, median=51.8421, max=51.8421
#>   dist   : available for 0/1 target(s)
```

For NP fits, `weighted_mean_y` equals the local mean estimate when the
weights are normalized.

## CATE Overlap

``` r

demo_cate <- cate_demo_data(n = 150, p = 2, design = "hetero", seed = 2)

fit_cate <- cate_dnn(
  data = demo_cate$data,
  x_cols = c("X1", "X2"),
  y_col = "Y",
  t_col = "D",
  m0_col = "m0_oracle", m1_col = "m1_oracle",
  e_col = "e_oracle",
  x0 = demo_cate$x0,
  method = "exact",
  k = 1L,
  s = 4L,
  se = FALSE,
  return = c("W", "dx", "psi", "nuisance")
)

d_cate <- diagnostics(fit_cate)
d_cate[, c(
  "target_id",
  "n_eff",
  "treated_mass",
  "control_mass",
  "treated_share",
  "weighted_mean_propensity",
  "treated_ess", "control_ess", "weighted_propensity_q10",
  "weighted_propensity_q90", "top5_score_share"
)]
#> <dnn_diagnostics>
#>   targets: 1
#>   n_eff  : min=64.9674, median=64.9674, max=64.9674
#>   dist   : available for 0/1 target(s)
#>   treated share: 0.429 to 0.429
#>   treated ESS: 28.4 to 28.4
#>   control ESS: 36.6 to 36.6
#>   propensity q10: 0.37 to 0.37
#>   propensity q90: 0.607 to 0.607
#>   top-five absolute score-contribution share: 0.217 to 0.217
```

Overlap diagnostics are available when the fit was constructed from
treatment and propensity inputs. If
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
was called with only a direct `pseudo_outcome`, overlap fields remain
present but are marked unavailable.

Printing summarizes available arm information, propensity quantiles and
top-five absolute score-contribution shares; the data frame retains each
profile’s values. CLATE distinguishes instrument support, clipped
propensity mass and squared jackknife-contribution concentration. These
quantities have different definitions and should not be interpreted
interchangeably. Clipping prevents undefined arithmetic but does not
establish overlap. The reported order fraction is descriptive and
triggers no new warning policy.

## Experimental PLQR conditioning

All PLQR functionality, including these diagnostics and its specialized
plots, is experimental and outside the CATE and regression focus of this
release.

``` r

demo_plqr <- plqr_demo_data(n = 160, p = 2, p_d = 2, tau = 0.5, seed = 3)

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
  x0 = demo_plqr$x0,
  method = "exact",
  k = 1L,
  s = 4L,
  se = TRUE,
  return = c("W", "dx", "score", "ess")
)

d_plqr <- diagnostics(fit_plqr)
d_plqr[, c(
  "target_id",
  "n_eff",
  "gram_min_eig",
  "gram_condition_number",
  "jacobian_min_singular_value",
  "jacobian_condition_number",
  "theta_vcov_available"
)]
```

The PLQR Gram diagnostic is based on the localized residual-regressor
matrix `D - m_D(X)`. It is an identification and conditioning
diagnostic, not a local bread estimate.

## Tuning Metadata

When a fit uses `tuning = ...`, the diagnostics table exposes compact
tuning metadata without duplicating the full candidate grid into every
row.

``` r

fit_tuned <- np_dnn(
  data = demo_np$data,
  x_cols = c("X1", "X2"),
  y_col = "Y",
  x0 = demo_np$x0,
  method = "exact",
  k = 1L,
  tuning = tuning_dnn_paper(lambda = 1.25),
  se = FALSE,
  return = c("W", "tuning")
)

d_tuned <- diagnostics(fit_tuned)
d_tuned[, c(
  "target_id",
  "tuning_name",
  "tuning_default_s",
  "tuning_grid_n_candidates",
  "tuning_lambda_default"
)]
#> <dnn_diagnostics>
#>   targets: 1
#>   dist   : available for 0/1 target(s)

attr(d_tuned, "tuning_grid")
#>   candidate_id computed_key is_default feasible reason status lambda  s
#> 1         s_45         s_45       TRUE     TRUE     ok     ok     NA 45
#> 2   lambda_0p3         s_10      FALSE     TRUE     ok     ok   0.30 10
#> 3   lambda_0p4         s_14      FALSE     TRUE     ok     ok   0.40 14
#> 4   lambda_0p5         s_18      FALSE     TRUE     ok     ok   0.50 18
#> 5  lambda_0p75         s_27      FALSE     TRUE     ok     ok   0.75 27
```

Use `attr(d_tuned, "tuning_grid")` when you need the full candidate
table.

## Availability Rules

Diagnostics avoid hidden recomputation.

- If `W` is missing,
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
  errors and asks for a refit with weights.
- If `dx` is missing, distance fields are returned as `NA` with
  `distance_diagnostics_available = FALSE`.
- If CATE treatment/propensity state is unavailable, overlap fields are
  present but marked unavailable.
- If PLQR residual-regressor state is unavailable, Gram fields are
  present but marked unavailable.

This behavior keeps diagnostics predictable and prevents the diagnostics
layer from silently rerunning expensive estimator construction.

## Related Articles

- [Notation and package
  conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md):
  [`vignette("notation-and-conventions", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md)
  distinguishes similarly named fields and moment- from
  coefficient-scale uncertainty.
- [Visualising
  fits](https://JakobJuergens.github.io/RdistNN/articles/visualising-fits.md):
  [`vignette("visualising-fits", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/visualising-fits.md)
  turns diagnostic fields into reusable plot panels.
- [Conditional
  Z-estimation](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md):
  [`vignette("conditional-z-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md)
  supplies the estimator workflows that create these objects.
