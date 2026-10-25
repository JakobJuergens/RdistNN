# Diagnostics and Reporting

[`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
describes localization and estimator-specific support from a fitted
object without refitting nuisance models. The examples below show the
returned tables and how to interpret them. Diagnostics do not verify the
paper’s inference conditions; distance-tie averaging can change support
and concentration, and inference with persistent distance ties is not
established.

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
as.data.frame(d_np)[, c(
  "target_id",
  "n_eff",
  "weighted_mean_y",
  "weighted_sd_y",
  "weighted_y_q50"
)]
#>          target_id    n_eff weighted_mean_y weighted_sd_y weighted_y_q50
#> target_1  target_1 51.84209         0.08698       0.54151      0.1091348
```

For NP fits, `weighted_mean_y` equals the local mean estimate when the
weights are normalized.

| Question | Fields to inspect | Interpretation |
|----|----|----|
| Is weight concentrated? | `n_eff`, `max_weight`, top-weight shares | Concentration, not a count of neighbors; signed TDNN concentration is not probability-weight ESS. |
| How far does localization reach? | Rank radii and weighted distance summaries | Depend on the selected metric; CATE/CLATE physical radii use original Euclidean coordinates. |
| Are both binary arms represented? | Arm mass and within-arm ESS | CATE uses treatment arms; CLATE uses instrument arms. Read mass and concentration together. |
| How dispersed is the local outcome? | `weighted_sd_y` and outcome quantiles | Outcome spread, not uncertainty about the estimated mean. Negative TDNN weights do not define these probability summaries. |
| Is ratio inference available? | CLATE components, first stage and statuses | Inspect these alongside the ratio; see the [paired-inference guide](https://JakobJuergens.github.io/RdistNN/articles/conditional-late-inference.md). |

For NP/TDNN and CATE, [`nobs()`](https://rdrr.io/r/stats/nobs.html)
counts training observations. Candidate summaries use minimum ESS and
maximum weight/concentration over targets; target rows retain each
profile’s values.

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
as.data.frame(d_cate)[, c(
  "target_id",
  "n_eff",
  "treated_mass",
  "control_mass",
  "treated_share",
  "weighted_mean_propensity",
  "treated_ess", "control_ess", "weighted_propensity_q10",
  "weighted_propensity_q90", "top5_score_share"
)]
#>          target_id    n_eff treated_mass control_mass treated_share
#> target_1  target_1 64.96742    0.4291109    0.5708891     0.4291109
#>          weighted_mean_propensity treated_ess control_ess
#> target_1                0.4935985    28.38508    36.59798
#>          weighted_propensity_q10 weighted_propensity_q90 top5_score_share
#> target_1               0.3696287               0.6065485        0.2167448
```

Overlap diagnostics are available when the fit was constructed from
treatment and propensity inputs. If
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
was called with only a direct `pseudo_outcome`, overlap fields remain
present but are marked unavailable.

`treated_ess` and `control_ess` normalize weights within their treatment
arm before measuring concentration. An arm with one positive-weight
observation has ESS one even if that observation receives tiny total
weight; an arm with zero weight has `NA` ESS. Read these fields
alongside `treated_mass` and `control_mass`: arm ESS describes
concentration within an arm, while arm mass is the total localization
weight assigned to that arm.

Printing summarizes available arm information, propensity quantiles and
top-five absolute score-contribution shares; the data frame retains each
profile’s values. CLATE distinguishes instrument support, clipped
propensity mass and squared jackknife-contribution concentration. These
quantities have different definitions and should not be interpreted
interchangeably. Clipping prevents undefined arithmetic but does not
establish overlap. The reported order fraction is descriptive, rather
than a calibrated validity threshold.

## Experimental PLQR conditioning

PLQR diagnostics and specialized plots remain experimental.

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
as.data.frame(d_plqr)[, c(
  "target_id",
  "n_eff",
  "gram_min_eig",
  "gram_condition_number",
  "jacobian_min_singular_value",
  "jacobian_condition_number",
  "theta_vcov_available"
)]
#>          target_id    n_eff gram_min_eig gram_condition_number
#> target_1  target_1 69.34251     0.286607              1.333248
#>          jacobian_min_singular_value jacobian_condition_number
#> target_1                   0.5965432                  1.274692
#>          theta_vcov_available
#> target_1                 TRUE
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
as.data.frame(d_tuned)[, c(
  "target_id",
  "tuning_name",
  "tuning_default_s",
  "tuning_grid_n_candidates",
  "tuning_lambda_default"
)]
#>          target_id       tuning_name tuning_default_s tuning_grid_n_candidates
#> target_1  target_1 paper_undersmooth               45                        5
#>          tuning_lambda_default
#> target_1                  1.25

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

- Missing `W` prevents NP/TDNN, CATE and PLQR target diagnostics. Stored
  grid summaries and paired-ratio/CLATE diagnostics have their own
  retained state.
- Rank radii and weighted metric distances require retained `dx`;
  otherwise those fields are `NA` with
  `distance_diagnostics_available = FALSE`. CATE/CLATE physical radii
  can still use retained observations in original units.
- If CATE treatment/propensity state is unavailable, overlap fields are
  present but marked unavailable.
- If PLQR residual-regressor state is unavailable, Gram fields are
  present but marked unavailable.

The
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
states the conditions for inference. These diagnostics do not check all
of them.

## Using diagnostics in scripts

Warnings summarize one issue family per public call. Independent reports
or plots can warn again. Base R warning controls and calling handlers
apply, including `options(warn = 2)`. To silence only the tie advisory,
use
`suppressWarnings(fit_call, classes = "rdistnn_distance_ties_warning")`.

Paired-ratio/CLATE calls retain every target, including unavailable
ratios or intervals, and point to `status` and `components`. Negative
denominators alone are descriptive; deliberately disabled uncertainty
does not warn. Ordinary mean fits instead stop on nonfinite computed
estimates and request rescaling. Dashboards may omit absent optional
extras; invalid options, inconsistent state and unexpected rendering
failures stop.

| Cause | Condition identity |
|----|----|
| Invalid input | `rdistnn_input_error` |
| Invalid extension return value | `rdistnn_callback_error` |
| Computed numerical-range failure | `rdistnn_numerical_error` |
| Missing optional extra | `rdistnn_extra_unavailable` |
| Invalid retained state, requiring refitting | `rdistnn_state_error` |
| Package warning | `rdistnn_warning`, with a family such as `rdistnn_distance_ties_warning` or `rdistnn_partial_result_warning` |

Owned conditions carry a semantic `code` and applicable public
`argument`, expected/actual facts and bounded `operation` name, without
storing the full call or its data. Handle recognized identities rather
than parsing English messages. Unexpected callback, native and resource
errors retain their identities. Unknown character `return` entries warn
and are ignored; missing or noncharacter selectors are input errors.

### Read notice facts without displaying warnings

`fit$meta$target_notices` retains `target_id`, `distance_ties`,
`weight_concentration` and `weak_localization`; `meta$n_tied_targets`
and `meta$n_weak_targets` retain totals. These fields are available even
with `return = character()` and survive serialization. Candidates
describe their selected weights. Concentration is `1 / sum(w^2)` from
final observation weights; for signed TDNN it is descriptive rather than
probability-weight ESS. The weak-localization guide is at least `.9 * n`
with `n >= 2`, not a validity test. Complete current fitted state is
required before reuse.

## Related Articles

- [Notation and
  conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md):
  field meanings and uncertainty scales.
- [Visualising
  fits](https://JakobJuergens.github.io/RdistNN/articles/visualising-fits.md):
  reusable plot panels.
- [Experimental
  PLQR](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md):
  moment and coefficient uncertainty.
