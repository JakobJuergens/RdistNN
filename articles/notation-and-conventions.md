# Notation and Package Conventions

This is a guide to package arguments, array shapes and result fields.
Start with [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md)
for an analysis you can run. For the statistical setup, assumptions and
derivations, see the accompanying
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).

## Localization and DNN conventions

`X` is an `n`-by-`p` matrix of localization covariates. Use `x0` for a
length-`p` profile or `X0` for an `m`-by-`p` matrix with one profile per
row. Supply exactly one target argument. In a data-frame call, `x_cols`
selects and orders the coordinates; targets must follow that order.

Reporting coordinates can be a subset of the covariates used by an
external nuisance learner. Keep nuisance predictions aligned with the
training rows, even if the learner used additional covariates.

### Kernel order `s` and neighbor count `k`

- `s` is the subsample order. Larger `s` gives tighter localization.
- `k` counts neighbors inside each subsample; the default is `k = 1`.
- `s` does not count retained full-sample neighbors. The unaveraged
  exact rank table reaches rank `n - s + k` (`n - s + 1` at the default
  `k = 1`). A duplicate block crossing that edge can give later
  observations positive averaged weight.
- Jackknife uncertainty needs `s < n`; estimate-only fits can use
  `s = n`.

When omitted, `s` follows
[`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md)
with multiplier `lambda = 0.4`. [Kernel-order
tuning](https://JakobJuergens.github.io/RdistNN/articles/kernel-order-tuning.md)
explains explicit orders and candidate grids. [DNN
weights](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md)
covers rank and weight output. The paper gives the theoretical scope of
the default construction.

### Weights and effective sample size

`W` has observations in rows and targets in columns. For normalized
weights, `ess`, `neff` and diagnostic `n_eff` all mean `1 / sum(w^2)`.
This measures weight concentration; it is not a count of nonzero weights
or a distance radius. For NP/TDNN and CATE, `nobs(fit)` instead counts
training observations. TDNN can have signed weights, so positive-weight
overlap interpretations do not carry over to it.

### Numerical geometry and covariance

Fresh fits order exact represented distance ties between distinct
covariate vectors lexicographically by the original supplied columns.
Within each group of exactly identical localization covariate rows, each
observation receives the average of the group’s original rank weights.
The block’s total mass is preserved; individual weights, support and ESS
can change. No distance tolerance groups nearby values, and outcomes,
treatments and nuisances do not define groups.

Each jackknife deletion removes the actual observation, updates ranks,
and averages the `n-1` weight table over every remaining block. This
makes estimates and covariance invariant to aligned row permutations up
to numerical precision. The package warns when duplicate groups occur
because the modification is not covered by the paper’s theory. Refit
tied saved objects created before this policy; rank-only utilities
without covariates do not infer duplicate groups.

Saved prediction reuses the fitted training geometry. Fast jackknife
calculations also freeze it; NP/TDNN/CATE `jackknife = "loo"` refits
sample-dependent geometry for each deletion. The two calculations can
differ. See
[`?distance_mahalanobis`](https://JakobJuergens.github.io/RdistNN/reference/distance_mahalanobis.md)
for the fitted metric and rank-deficient cases. Numerical engine details
belong to the maintainer reference linked from
[CONTRIBUTING](https://github.com/JakobJuergens/RdistNN/blob/master/CONTRIBUTING.md).

## Matrix and data-frame interfaces

| Quantity | Data-frame input | Matrix/vector input |
|----|----|----|
| Localization covariates | `data`, `x_cols` | `X` |
| Reporting profiles | `x0` or `X0` | `x0` or `X0` |
| Regression outcome | `y_col` or `y` | `Y` |
| CATE signal | nuisance column names or `pseudo_outcome` | `pseudo_outcome` |
| PLQR regressors | `D_cols` | `D` |
| PLQR nuisances | `g_col`, `mD_cols`, `density_col` | `g`, `mD`, `density` |

Tibbles use the data-frame route. Column selection is explicit; the
package does not guess variable names. These equivalent calls use a data
frame with numeric `X1`, `X2` and `Y` columns:

``` r

df <- np_demo_data(n = 60, p = 2, seed = 1)$data
fit_table <- np_dnn(data = df, x_cols = c("X1", "X2"), y_col = "Y", x0 = c(0, 0))
fit_matrix <- np_dnn(X = as.matrix(df[c("X1", "X2")]), Y = df$Y, x0 = c(0, 0))
all.equal(coef(fit_table), coef(fit_matrix))
#> [1] TRUE
```

## Estimator targets and signals

| Function | Supplied input beyond covariates/targets | Main result |
|----|----|----|
| [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md) | Outcome or scalar signal | `mu`: conditional mean |
| [`np_tdnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_tdnn.md) | Outcome and two orders `s1 < s2` | `mu`: signed two-scale regression estimate |
| [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md) | AIPW nuisance predictions or supplied pseudo-outcome | `tau`: conditional treatment effect |
| [`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md) | Outcome, binary instrument/treatment and five nuisance predictions | `theta`: conditional complier effect; paired `components` |
| [`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md) | Numerator and denominator signals | `theta`: ratio of local averages |
| [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md) | Regressors, nuisance predictions and `theta_eval` | `moment` at the supplied coefficient |
| [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md) | A moment model and its data/nuisances | Solved `theta` and convergence diagnostics |

A **score** is an observation-level contribution. A **moment** is its
local weighted average. The experimental PLQR wrapper evaluates a
moment; it does not solve for a coefficient. Experimental GMM does solve
a moment problem.

## CATE nuisance conventions

For AIPW, supply `y_col`, `t_col`, `m0_col`, `m1_col` and `e_col`:
outcome, binary treatment, control mean, treated mean and treatment
propensity. Alternatively supply `pseudo_outcome`, which is used
unchanged and is the CATE route for matrix inputs. Nuisance models are
fitted outside the package.

Propensities must be finite and in `[0, 1]`. The AIPW route clips them
to `[clip, 1 - clip]`, with default `clip = 1e-8`; clipping does not
establish overlap. A supplied signal alone cannot supply
treatment-overlap diagnostics. Earlier pooled-mean `m_col` inputs have
been removed; both arm means are required for AIPW. See the [migration
guidance](https://JakobJuergens.github.io/RdistNN/articles/getting-started.html#migrating-earlier-cate-inputs).

## CLATE conventions

`A` is the binary instrument and `D` is treatment. The five predictions
`e`, `m_y0`, `m_y1`, `m_d0`, `m_d1` are the instrument propensity and
outcome/treatment means in each **instrument arm**. `one_sided = TRUE`
selects the structural restriction that treatment is unavailable without
assignment; it permits omitting `m_d0`. Observing `D <= A` does not by
itself justify that restriction.

[`iv_augmented_signals()`](https://JakobJuergens.github.io/RdistNN/reference/iv_augmented_signals.md)
constructs paired outcome and treatment signals. The fit divides their
common-weight local averages. Augmented treatment signals may be
negative and are not compliance indicators. See the [CLATE
guide](https://JakobJuergens.github.io/RdistNN/articles/conditional-late.md)
for complete calls and [paired
inference](https://JakobJuergens.github.io/RdistNN/articles/conditional-late-inference.md)
for covariance and status outputs. Ordinary residualized IV need not
have the same causal target.

## Experimental PLQR target and nuisance conventions

All PLQR interfaces remain experimental. `theta_eval` is supplied by the
user; `theta` repeats it. The principal inputs are:

| Input | Meaning |
|----|----|
| `D` | Partially linear regressors |
| `g` | Nonparametric shape nuisance; not the entire conditional quantile |
| `mD` | Conditional means used to residualize the regressors |
| `density` | Conditional outcome density at the model threshold |

The package does not fit those nuisances or solve a local coefficient
root. See [Experimental
PLQR](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md)
for calls and
[`?plqr_dnn`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
for the full contract.

## Moment diagnostics versus coefficient inference

PLQR uncertainty has two distinct scales:

| Output                     | Interpretation                             |
|----------------------------|--------------------------------------------|
| `moment_se`, `moment_vcov` | Uncertainty of the localized moment        |
| `theta_se`, `theta_vcov`   | Density-backed coefficient sandwich output |
| `theta_omega`              | `theta_vcov` multiplied by `neff`          |

Coefficient intervals require `density`/`density_col`. Without it,
`moment_se` is a diagnostic and does not support
[`confint()`](https://rdrr.io/r/stats/confint.html) or prediction
intervals for coefficients. With density, `theta_meat = "analytic"` uses
score covariance; `theta_meat = "jackknife"` uses joint moment
covariance. These calculations do not turn a supplied coefficient into
an estimated root.

## Experimental GMM conventions

[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
uses a `gmm_model` to solve local moments and returns its own
coefficient covariance and convergence diagnostics. Its `J` is an
overidentification statistic; PLQR’s `J`/`jacobian` is a derivative
matrix. See [Local IV
estimation](https://JakobJuergens.github.io/RdistNN/articles/local-iv-estimation.md)
for the experimental GMM workflow.

## Reading fitted objects

Stored fields and available methods differ by estimator:

| Class | Estimate fields | Uncertainty |
|----|----|----|
| `np_result` / `tdnn_np_result` | `mu` | `se`; joint target covariance from [`vcov()`](https://rdrr.io/r/stats/vcov.html) |
| `cate_result` | `tau` | `se`; joint target covariance from [`vcov()`](https://rdrr.io/r/stats/vcov.html) |
| `affine_ratio_result` / `clate_result` | `theta`, `components` | `se`, `component_vcov`, `theta_vcov` |
| `plqr_result` | `moment`, supplied `theta_eval` | Moment- and coefficient-scale fields above |
| `gmm_result` | Solved `theta` | `se`, `vcov` |

CATE and NP/TDNN covariance uses common deleted observations across
targets and requires an active fitted jackknife. It is calculated from
retained state, including serialized fits without optional extras. Both
axes follow [`coef()`](https://rdrr.io/r/stats/coef.html) order.
Prespecified contrasts still need the corresponding joint inference
assumptions.

CLATE covariance lists are marginal across targets.
`component_vcov[[j]]` is the paired two-by-two component covariance and
`theta_vcov[[j]]` is the ratio variance. These are already on the
estimator scale: do not divide by `neff` again. `theta_omega[[j]]`
instead equals `neff[j] * theta_vcov[[j]]`.

`return` extras such as `W`, `dx`, `ess` and `tuning` add retained
information without changing the main estimate. Use
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
as documented for the fitted class. Missing uncertainty remains missing.
Finite ratio estimates can have unavailable Wald intervals; inspect
`status` and `wald_available` before reporting them.

Saved objects retain observation-level numeric inputs, outcomes or
signals; apply the underlying data’s sharing restrictions. [Diagnostics
and
reporting](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md)
and
[visualisation](https://JakobJuergens.github.io/RdistNN/articles/visualising-fits.md)
cover inspection and presentation. Intervals are pointwise; none of
these output conventions establishes the paper’s assumptions or
simultaneous confidence bands.
