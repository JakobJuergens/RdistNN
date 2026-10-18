# Notation and Package Conventions

## Purpose

`RdistNN` uses the notation of its companion papers where that notation
has a direct computational meaning. This article collects the
translation in one place. It is intended to answer three questions that
otherwise recur across help pages:

1.  Which mathematical object does a package argument represent?
2.  Is the object supplied by the user, computed by the package, or
    merely retained as a diagnostic?
3.  Does an uncertainty field describe a localized moment or an
    estimated target coefficient?

The notation is shared across several estimators, but the estimators do
not all perform the same operation. In particular,
[`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md)
and
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
localize scalar signals, the experimental
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
evaluates a vector moment at a supplied coefficient, and the
experimental
[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
front end solves a local moment problem.

This is the collection’s crosswalk and navigation article. It explains
notation and conventions that cut across estimators; the linked articles
remain the authoritative source for complete workflows, diagnostics,
extension patterns, and paper-specific replication.

## How the vignette collection fits together

The vignette collection has four complementary layers:

1.  Start with the package tour and the common DNN localization
    machinery.
2.  Move to an estimator article for the NP and CATE workflows or
    experimental PLQR and IV details.
3.  Consult the diagnostics and visualization articles when interpreting
    a fit.
4.  Use the extension and replication articles only when you need those
    more specialized interfaces.

This article sits between those layers. Return here when a symbol or
field has different meanings in two workflows, or when it is unclear
whether an object is an input, an intermediate quantity, a diagnostic,
or an inferential output.

## Find information by task

| If you want to … | Continue with … |
|----|----|
| Take a first tour of weights and the main estimators | [Getting started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md): [`vignette("getting-started", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md) |
| See NP, CATE, and experimental PLQR calls | [Conditional Z-estimation](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md): [`vignette("conditional-z-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md) |
| Estimate a binary-instrument complier effect from supplied predictions | [Conditional LATE](https://JakobJuergens.github.io/RdistNN/articles/conditional-late.md) and [paired inference](https://JakobJuergens.github.io/RdistNN/articles/conditional-late-inference.md) |
| Understand exact and asymptotic DNN weights, support, and ranks | [DNN weights](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md): [`vignette("dnn-weights", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md) |
| Choose `s` with the paper rule, registry, or tuning grid | [Kernel-order tuning](https://JakobJuergens.github.io/RdistNN/articles/kernel-order-tuning.md): [`vignette("kernel-order-tuning", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/kernel-order-tuning.md) |
| Interpret diagnostic names and their availability | [Diagnostics and reporting](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md): [`vignette("diagnostics-reporting", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md) |
| Plot fitted targets, weight geometry, and overlap | [Visualising fits](https://JakobJuergens.github.io/RdistNN/articles/visualising-fits.md): [`vignette("visualising-fits", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/visualising-fits.md) |
| Define a custom distance, rank-weight map, or score | [Extending RdistNN](https://JakobJuergens.github.io/RdistNN/articles/extending.md): [`vignette("extending", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/extending.md) |
| Work directly with generalized U-statistic jackknife tools | [Jackknife generalized U-statistics](https://JakobJuergens.github.io/RdistNN/articles/jackknife-gen-ustats.md): [`vignette("jackknife-gen-ustats", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/jackknife-gen-ustats.md) |
| Use the experimental local IV-GMM workflow | [Local IV estimation](https://JakobJuergens.github.io/RdistNN/articles/local-iv-estimation.md): [`vignette("local-iv-estimation", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/local-iv-estimation.md) |
| Reproduce the current affine simulations and SIPP application | [Affine paper replication](https://JakobJuergens.github.io/RdistNN/articles/paper-affine-replication.md): [`vignette("paper-affine-replication", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/paper-affine-replication.md) |

## Quick package-to-paper map

| Package name | Paper notation | Role |
|----|----|----|
| `X` | $`X_i`$ | Localization covariates for observation $`i`$ |
| `x0` | $`x`$ or $`x_0`$ | One target point |
| `X0` | $`(x_1,\ldots,x_M)^\top`$ | Matrix of target points, one per row |
| `s` | $`s_n`$ | DNN kernel order, or subsample size in the complete-statistic representation |
| `k` | $`k`$ | Fixed nearest-neighbor count inside each size-`s` subsample |
| `W` | $`[w_i(x_j)]`$ | Observation-by-target DNN weight matrix |
| `n_eff` / `neff` / `ess` | $`n_{\mathrm{eff}}(x)`$ | Effective local sample size $`1/\sum_i w_i(x)^2`$ |
| `lambda` | $`\lambda_s`$ | Finite-sample multiplier in the paper tuning rule for `s` |
| `theta_eval` | $`\theta`$ at which $`m(x;\theta,\eta)`$ is evaluated | Supplied PLQR coefficient; not solved for by [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md) |
| `theta` | $`\theta(x)`$ | Target coefficient; supplied in PLQR and solved in experimental GMM |
| `g` | $`g_0(U)`$ | PLQR nonparametric shape nuisance |
| `mD` | $`m_D(U)=\mathbb E[D\mid U]`$ | PLQR residualizing nuisance |
| `density` | $`q_{0,x}(D,Z)`$ | Threshold-density input for the PLQR Jacobian |
| `score` | $`\psi_i`$ | Observation-level contribution before localization |
| `moment` | $`\sum_i w_i(x)\psi_i`$ | Localized empirical moment |

Names vary slightly when R object conventions require it. For example,
`plqr_result` uses `neff`, whereas a requested common weight diagnostic
is called `ess`. These fields represent the same inverse-Herfindahl
effective sample size when the weight columns are normalized.

## Localization and DNN conventions

### Training covariates and targets

`X` is the $`n\times p`$ matrix used to rank observations around a
target. The package uses `x0` for one length-$`p`$ target and `X0` for
an $`M\times p`$ matrix of targets. Supply exactly one of `x0` and `X0`.
Internally, both are resolved to a target matrix whose rows index
evaluation points.

In the data-frame interface, `x_cols` selects the columns that form `X`.
The values in `x0` or rows of `X0` must follow that same column order.
When column names are available, fitted objects retain them in the
returned `x0` matrix.

The localization state $`X`$ need not equal the full adjustment state
used by a nuisance function. The conditional Z-estimation paper often
writes $`U=(X,Z)`$. The package localizes on `X`, while supplied
nuisances such as PLQR `g` and `mD` may have been estimated using the
larger state $`U`$.

### Kernel order `s` and neighbor count `k`

The two integer controls have different jobs:

- `s` is the DNN kernel order. In the complete-statistic interpretation,
  it is the subsample size over which the nearest-neighbor kernel is
  averaged.
- `k` is the number of nearest observations selected *inside each
  size-`s` subsample*. The companion paper’s stated theory uses the
  default `k = 1L`.

Consequently, `s` is not an ordinary request to keep the nearest `s`
observations. Exact DNN weights can have support beyond rank `s`; for
fixed `k`, their finite-sample support reaches rank `n - s + k`. The
optional `neighbors` return records the first `k` observations used
inside a subsample, not the complete support of the averaged DNN weight
vector.

When `s` is omitted,
[`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md)
supplies the paper rule

``` math
s = \left\lfloor
\lambda_s n^{(k_x+1)/(k_x+2)}
\right\rfloor,
```

with package argument `lambda`, default value `0.4`, and $`k_x=\dim(X)`$
unless `kx_override` is set. The argument is called `lambda` rather than
`kappa` so it does not collide with unrelated uses of $`\kappa`$ in the
theory.

### Weights and effective sample size

For one target, let $`w_i(x)`$ denote the normalized DNN weight on
observation $`i`$. For multiple targets, `W[i, j]` is $`w_i(x_j)`$. With
normalized weights,

``` math
n_{\mathrm{eff}}(x)
=
\left\{\sum_{i=1}^n w_i(x)^2\right\}^{-1}.
```

The same quantity appears under three public names:

- `ess` when requested as a common weight extra;
- `neff` in PLQR and GMM result objects;
- `n_eff` in tidy diagnostics and paper-facing reporting.

It is a weight-concentration measure, not a count of nonzero weights.
Physical radius diagnostics answer a different question: how far from
`x0` one must move to collect a specified fraction of weight mass.

Continue with [DNN
weights](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md)
for the weight algorithms and returned rank diagnostics, or
[kernel-order
tuning](https://JakobJuergens.github.io/RdistNN/articles/kernel-order-tuning.md)
for the tuning rule and candidate-grid workflow.

## Matrix and data-frame interfaces

High-level estimators support two outer interfaces. The data-frame
interface is usually easier to read; the matrix interface is useful in
programmatic code. After validation, both resolve to the same numeric
matrices and vectors.

| Concept | Data-frame interface | Matrix interface |
|----|----|----|
| Localization covariates | `data` plus `x_cols` | `X` |
| One or many targets | `x0` or `X0` | `x0` or `X0` |
| NP outcome | `y_col` or `y` | `Y` |
| CATE signal | nuisance column names or `pseudo_outcome` | `pseudo_outcome` |
| PLQR regressors | `D_cols` | `D` |
| PLQR nuisances | `g_col`, `mD_cols`, `density_col` | `g`, `mD`, `density` |

Tibbles enter through the data-frame route; there is no separate tibble
API. The package does not guess columns through nonstandard evaluation.
Column names and their order are explicit parts of the call.

``` r

# Data-frame grammar
np_dnn(data = df, x_cols = c("X1", "X2"), y_col = "Y", x0 = c(0, 0))

# Equivalent matrix grammar
np_dnn(X = as.matrix(df[c("X1", "X2")]), Y = df$Y, x0 = c(0, 0))
```

Continue with [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md)
for a compact tour of the preferred interface, and [Conditional
Z-estimation](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md)
for side-by-side estimator calls.

## Estimator targets and signals

The high-level wrappers localize with DNN weights, but they consume and
return different statistical objects.

| Wrapper | Localized input | Main returned target | Interpretation |
|----|----|----|----|
| [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md) | Outcome or scalar signal $`Y_i`$ | `mu` | Local conditional mean |
| [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md) | CATE pseudo-outcome $`\psi_i`$ | `tau` | Local conditional average treatment effect |
| [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md) | Vector PLQR score $`\psi_i(\theta_{\mathrm{eval}})`$ | `moment`; supplied `theta` | Moment evaluation at a user-supplied coefficient |
| [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md) | Model-defined vector moments $`g_i(\theta)`$ | solved `theta` | Experimental local GMM estimate |
| [`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md) | Paired augmented reduced-form and first-stage signals | `theta`, `components` | Experimental conditional complier effect under the IV assumptions |
| [`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md) | Two supplied scalar signals | `theta`, `components` | Experimental ratio of common-weight local averages |

The word **score** means an observation-level contribution. The word
**moment** means its DNN-weighted average at a target. For scalar NP and
CATE workflows, the weighted score is already the target estimate. For
vector conditional moments, a localized moment must additionally be
connected to a coefficient through evaluation or root solving.

## Where package objects come from

Separating supplied quantities from computed quantities prevents several
common misreadings of a fitted object. The table describes the
high-level wrappers; lower-level weight and score tools expose more of
the intermediate machinery.

| Wrapper | Supplied by the user | Computed by `RdistNN` | Main retained output |
|----|----|----|----|
| [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md) | `X`, targets, and an outcome or signal | DNN weights and their weighted mean | `mu`, and `se` when requested |
| [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md) | `X`, targets, and either nuisance evaluations or a pseudo-outcome | CATE signal when needed, then its DNN localization | `tau`, and `se` when requested |
| [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md) | `X`, targets, `D`, `theta_eval`, `g`, `mD`, and optionally `density` | Scores and their local moment; Jacobian and sandwich pieces when inputs permit | `moment`; supplied `theta`; separate moment- and coefficient-scale uncertainty fields |
| [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md) | `X`, targets, a `gmm_model`, observed data, and model-specific nuisances | Local moments, two-step root solving, weighting matrix, and sandwich pieces | solved `theta`, `se`, `vcov`, and overidentification diagnostics |
| [`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md) | `X`, targets, `Y,D,A` and five nuisance predictions | Both augmented signals, paired local averages and frozen covariance | `components`, `theta`, `se`, per-target `status` |
| [`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md) | `X`, targets, numerator and denominator signals | Common-weight averages, ratio and paired frozen covariance | `components`, `theta`, `se`, per-target `status` |

The wrappers do not fit CATE, CLATE, PLQR, or built-in IV-GMM nuisances.
Those functions receive already evaluated nuisance vectors or matrices.
Also note the asymmetry in `theta`:
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
repeats a supplied evaluation value, whereas
[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
estimates a local root.

Continue with [Conditional
Z-estimation](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md)
for the estimator workflows, or [Extending
RdistNN](https://JakobJuergens.github.io/RdistNN/articles/extending.md)
for the lower-level distance, rank-weight, and score layers.

## Names that change meaning across estimators

Some short names are conventional within a model but are not
package-wide synonyms. Interpret them together with the fitted-object
class.

| Name | Meaning in one context | Meaning elsewhere or important distinction |
|----|----|----|
| `theta` | Supplied PLQR evaluation coefficient, copied into a `plqr_result` | A solved local GMM coefficient, or a paired local ratio in an affine-ratio/CLATE result |
| `g` | The PLQR nonparametric shape nuisance $`g_0(U)`$ | `g_obs` is the observation-by-moment contribution matrix in GMM, not a nuisance named `g` |
| `J` / `jacobian` | `jacobian` is the derivative or bread matrix used for coefficient inference | `J` in `gmm_result` is the scalar local overidentification statistic |
| `se` | Target-scale standard error for NP, CATE, and GMM results | PLQR uses `moment_se` and `theta_se` so the inferential scale is explicit |
| `ess` / `neff` / `n_eff` | Three public spellings of inverse weight concentration | None is a physical neighbor count or a radius diagnostic |
| `score` / `moment` | An observation-level contribution versus its localized weighted average | For scalar mean/CATE estimation the average is the target; for vector Z-estimation it may still require evaluation or solving |

## CATE nuisance conventions

[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
accepts one of three signal specifications.

### Full AIPW signal

With outcome nuisances $`m_d(u)=\mathbb E[Y\mid D=d,U=u]`$ and
propensity $`e(u)=\Pr(D=1\mid U=u)`$, the package constructs

``` math
\psi_i
=
m_1(U_i)-m_0(U_i)
+\frac{D_i\{Y_i-m_1(U_i)\}}{e(U_i)}
-\frac{(1-D_i)\{Y_i-m_0(U_i)\}}{1-e(U_i)}.
```

The data-frame arguments are `m0_col`, `m1_col`, and `e_col`, together
with `y_col` and `t_col`.

### Robinson-form signal

With $`m(u)=\mathbb E[Y\mid U=u]`$, the alternative signal uses `m_col`
and `e_col`:

``` math
\psi_i
=
\frac{\{Y_i-m(U_i)\}\{D_i-e(U_i)\}}
     {e(U_i)\{1-e(U_i)\}}.
```

### Precomputed pseudo-outcome

Users may instead provide `pseudo_outcome` directly. This is the only
CATE signal route in the matrix interface. It is also the cleanest route
when first-stage estimation is managed by another package or a
cross-fitting pipeline.

The nuisance functions are supplied externally;
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
does not fit outcome or propensity models. It validates and clips
supplied propensity values using `clip`, constructs the requested
signal, and localizes that signal. Orthogonality of the score does not
by itself certify the stability of an arbitrary same-sample nuisance
learner.

Continue with [Conditional
Z-estimation](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md)
for a worked CATE fit. For the paper-specific data-generating process
and Monte Carlo harness, use [Affine paper
replication](https://JakobJuergens.github.io/RdistNN/articles/paper-affine-replication.md).

## Experimental CLATE conventions

CLATE uses `A` for the binary instrument and `D` for binary treatment.
Reporting uses `X`; nuisance predictions adjust for all pre-instrument
variables in $`V=(X,Z)`$. `Z` denotes additional adjustment covariates
in this notation. The five predictions `e,m_y0,m_y1,m_d0,m_d1` are
instrument propensity and outcome/treatment means indexed by
**instrument arm**. `one_sided = TRUE` explicitly selects $`D(0)=0`$,
permitting `m_d0 = NULL`. The observed check `D <= A` does not establish
that structural restriction.

[`iv_augmented_signals()`](https://JakobJuergens.github.io/RdistNN/reference/iv_augmented_signals.md)
returns `phi_y` and `phi_d`. In the affine score $`B\theta+b`$, CLATE
has $`B=\phi_D`$ and $`b=-\phi_Y`$. The fitted solution is the ratio of
their DNN-weighted averages. Individual augmented denominators may be
negative; they are not compliance indicators and are not clipped.
Ordinary residualized IV can introduce different stratum weights.

| Field | Scale or convention |
|----|----|
| `components` | Target-by-two matrix with `numerator`, `denominator` columns |
| `component_vcov[[j]]` | Paired 2-by-2 estimator covariance |
| `theta_vcov[[j]]` | 1-by-1 ratio variance, already on the estimator scale |
| `theta_omega[[j]]` | `neff[j] * theta_vcov[[j]]` |
| `se[j]` | Square root of the ratio variance |
| `status` | One row per reporting target, retained even for undefined ratios |

Covariance is marginal across targets. `ratio_status` distinguishes
`ok`, `zero_denominator`, `component_overflow` and `ratio_overflow`.
`negative_denominator` and `near_zero_denominator` are flags, with no
silent denominator replacement. `covariance_status` separately records
missing, positive, zero-score or numerically failed covariance
calculations. `wald_available` controls interval availability. An exact
zero denominator gives `NA` ratio/interval while retaining component
uncertainty; an exact zero score norm gives computational SE zero and
unavailable Wald bounds.

[`predict()`](https://rdrr.io/r/stats/predict.html) relocalizes the same
frozen signals at new reporting points.
[`ratio_confidence_set()`](https://JakobJuergens.github.io/RdistNN/reference/ratio_confidence_set.md)
instead accepts components and paired covariance and returns interval
unions, including rays, all-real and empty sets; unresolved numerical
geometry has a separate status. It does not replace
[`confint()`](https://rdrr.io/r/stats/confint.html) or establish uniform
weak-instrument coverage. The [CLATE inference
article](https://JakobJuergens.github.io/RdistNN/articles/conditional-late-inference.md)
gives the exact four-row example, scale conventions and statistical
conditions.

CLATE supports exact/asymptotic single-scale DNN with the same `method`,
`k`, `s`, `normalize` and `tuning_mode` controls as CATE and NP
regression. Two-scale TDNN is available for NP regression only.
Weighting choices do not by themselves verify the inference conditions
for a particular estimator.
[`summary()`](https://rdrr.io/r/base/summary.html) returns the common
coefficient table plus CLATE-specific `components` and `status` tables.
Its unavailable Wald bounds remain `NA`.

## Experimental PLQR target and nuisance conventions

All PLQR functionality is experimental, including its scores, helpers,
post-estimation methods, and plots. The initial release focuses on CATE
and nonparametric regression; the PLQR API and inference interfaces may
change.

For partially linear quantile regression, let $`D_i`$ be the partially
linear regressor and let $`U_i`$ denote the full adjustment state. At a
supplied coefficient `theta_eval`,
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
forms the observation-level score

``` math
\psi_i(\theta_{\mathrm{eval}})
=
\left[
\tau-
\mathbf 1\{Y_i\le D_i^\top\theta_{\mathrm{eval}}+g_0(U_i)\}
\right]
\left[D_i-m_D(U_i)\right].
```

The key names are deliberately literal:

- `theta_eval` is the coefficient at which the moment is evaluated.
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
  does not solve for a new coefficient. Its returned `theta` repeats
  `theta_eval` at each target and records
  `meta$theta_source = "supplied"`.
- `g` or `g_col` supplies $`g_0(U)`$, the nonparametric shape component.
  It is not itself a conditional quantile; the threshold is
  $`D^\top\theta+g_0(U)`$.
- `mD` or `mD_cols` supplies $`m_D(U)=\mathbb E[D\mid U]`$, used to form
  the residualized regressor $`\widetilde D=D-m_D(U)`$.
- `density` or `density_col` supplies the positive conditional density
  of $`Y`$ at the threshold. It enters the coefficient Jacobian, not the
  moment or its score-covariance meat.

These nuisances are supplied externally. They may depend on variables
outside `x_cols`; `RdistNN` receives their evaluated values and does not
fit or restrict their underlying models.

## Moment diagnostics versus coefficient inference

This distinction is the most important PLQR reporting convention.

A localized vector moment can have a componentwise sampling diagnostic
without supporting inference for a coefficient. In particular,
`moment_se` describes the components of `moment`; by itself it does not
estimate a coefficient standard error. Coefficient inference requires
the full sandwich

``` math
\widehat V_\theta(x)
=
\widehat J(x)^{-1}
\widehat S(x)
\widehat J(x)^{-\top}
/n_{\mathrm{eff}}(x).
```

The PLQR routes are:

| Requested route | Bread | Meat | Coefficient inference fields |
|----|----|----|----|
| No `density`, jackknife enabled | None | Componentwise moment diagnostic | `moment_se` only |
| `density`, `theta_meat = "analytic"` | Density-backed `jacobian` | Analytic `score_covariance` | `theta_se`, `theta_vcov`, `theta_omega` |
| `density`, `theta_meat = "jackknife"` | Density-backed `jacobian` | Full delete-1 `moment_vcov` | `moment_se`, `moment_vcov`, `theta_se`, `theta_vcov`, `theta_omega` |

For the analytic route,

``` math
\widehat S(x)
=
\tau(1-\tau)
\sum_i w_i(x)\widetilde D_i\widetilde D_i^\top,
\qquad
\widehat J(x)
=
\sum_i w_i(x)q_i\widetilde D_iD_i^\top.
```

For the jackknife-meat route, the full covariance `moment_vcov` replaces
the analytic meat while the same density-backed Jacobian remains. The
package does not reconstruct a coefficient covariance from the diagonal
entries in `moment_se`; doing so would discard cross-moment covariance.

The coefficient fields have distinct scales:

- `theta_vcov` is the finite-sample covariance used by
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  `predict(..., se.fit = TRUE)`;
- `theta_se` is `sqrt(diag(theta_vcov))` target by target;
- `theta_omega` is the effective-sample-size-scaled covariance,
  `theta_vcov * neff`.

Direct coefficient inference requires normalized weights because the
sandwich is defined for weight columns summing to one.

Continue with [Conditional
Z-estimation](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md)
for the full PLQR call and with [Diagnostics and
reporting](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md)
for Jacobian, conditioning, and field-availability diagnostics.

## Experimental GMM conventions

[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
is a separate, experimental front end. A `gmm_model` supplies a vector
moment function, coefficient names, starting logic, and—when
available—an analytic Jacobian. Unlike
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
solves for `theta` at each target using a two-step local GMM procedure.

Several names should not be transferred mechanically between the two
APIs:

- `g` is the PLQR shape nuisance, while `g_obs` in a `gmm_result` is the
  model-defined observation-level moment matrix.
- `jacobian` is the derivative of the local moment map. The scalar `J`
  in a `gmm_result` is instead the local overidentification statistic.
- `omega` is the estimated moment covariance entering the GMM weighting
  and sandwich calculation; `vcov` is the resulting coefficient
  covariance.
- `neff` remains the DNN effective sample size.

Nuisance values required by a supplied `gmm_model` are also supplied
externally. The package’s built-in doubly robust IV model validates
those columns but does not fit them.

Continue with [Local IV
estimation](https://JakobJuergens.github.io/RdistNN/articles/local-iv-estimation.md)
for the experimental IV-GMM workflow. Use [Extending
RdistNN](https://JakobJuergens.github.io/RdistNN/articles/extending.md)
when the task is to define reusable score machinery rather than fit the
built-in IV model.

## Reading fitted objects

The core result fields are intentionally estimator-specific.

| Class | Point or moment fields | Uncertainty fields | Useful metadata |
|----|----|----|----|
| `np_result` | `mu` | `se` | `x0`, `meta$k`, `meta$s`, `meta$method` |
| `cate_result` | `tau` | `se` | `x0`, optional `psi`, optional `nuisance` |
| `plqr_result` | `moment`, supplied `theta`, `theta_eval` | `moment_se` or sandwich fields described above | `meta$theta_source`, `meta$theta_meat_method`, `neff` |
| `gmm_result` | solved `theta`, `theta_first` | `se`, `vcov`; optional `omega`, `jacobian` | `J`, `df`, `pval_J`, `neff`, `convergence` |
| `affine_ratio_result` / `clate_result` | solved `theta`, two-column `components` | `component_vcov`, `theta_vcov`, `theta_omega`, `se` | `status`, `neff`, component roles and training declaration |

Optional fields requested through `return` are additive. For example,
asking for `W`, `dx`, `ess`, or `tuning` does not change the meaning of
the main estimate. Use
[`summary()`](https://rdrr.io/r/base/summary.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`predict()`](https://rdrr.io/r/stats/predict.html), and
[`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
where the fitted-object class documents those methods instead of
inferring meaning from similarly named internal fields.

Shared reporting functions dispatch on the estimation object. Use the
same calls below for CATE and CLATE; their score and uncertainty
calculations remain method-specific.

| Interface | CATE | CLATE |
|----|----|----|
| `summary(fit)$coefficients` | Estimate, SE, statistic, p-value and confidence bounds | Same columns; component/status tables additionally retained |
| [`plot_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_estimates.md), [`plot()`](https://rdrr.io/r/graphics/plot.default.html), [`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html) | Profile, forest or pointrange; zero reference | Same styles and arguments; unavailable ratio bounds stay missing |
| [`plot_weights()`](https://JakobJuergens.github.io/RdistNN/reference/plot_weights.md) | Profile, cumulative or kernel | Same views and retained-weight requirements |
| [`plot_overlap()`](https://JakobJuergens.github.io/RdistNN/reference/plot_overlap.md) | Treatment mass and propensity | Instrument mass and propensity |
| [`plot_influence()`](https://JakobJuergens.github.io/RdistNN/reference/plot_influence.md) | Weighted pseudo-outcome residual | Weighted paired score residual |
| [`plot_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/plot_tuning.md), [`candidate_fit()`](https://JakobJuergens.github.io/RdistNN/reference/candidate_fit.md), [`candidate_fits()`](https://JakobJuergens.github.io/RdistNN/reference/candidate_fit.md) | Frozen-state sensitivity and candidate views | Same workflow, retaining paired covariance |
| [`plot_dashboard()`](https://JakobJuergens.github.io/RdistNN/reference/plot_dashboard.md) | Common component panels | Same layout and optional-dependency policy |

The tidy
[`plot_data_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_data_estimates.md)
and
[`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
tables allow custom presentations without accessing plotting internals.
Distinct result classes retain their scientific meaning; CLATE is not a
`cate_result` subclass.

``` r

fit_cate$tau
fit_cate$se
diagnostics(fit_cate)

fit_plqr$moment
fit_plqr$theta_eval
fit_plqr$theta_se       # coefficient scale, when density was supplied
fit_plqr$moment_se      # moment scale; not a coefficient SE by itself
```

Continue with [Diagnostics and
reporting](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md)
to see which diagnostics are structural, requested, or
estimator-specific, and with [Visualising
fits](https://JakobJuergens.github.io/RdistNN/articles/visualising-fits.md)
to turn those fields into plots.

## Practical reading checklist

When translating between a paper display and an `RdistNN` call, check
the objects in this order:

1.  Identify the localization covariates `X` and target `x0` or `X0`.
2.  Record `s`, `k`, the weight method, and whether weights are
    normalized.
3.  Identify the observation-level signal or vector score being
    localized.
4.  Mark every nuisance as supplied or package-computed. Current CATE,
    CLATE, PLQR, and built-in GMM workflows expect nuisance values
    rather than fitting them.
5.  Determine whether the returned uncertainty is on the moment scale or
    the target-coefficient scale.
6.  Use `n_eff`/`neff`/`ess` as a scale and concentration diagnostic,
    not as a substitute for checking physical localization and
    first-stage assumptions.

The task table near the beginning of this article is the quickest route
onward. For most analyses, the next stop is [Conditional
Z-estimation](https://JakobJuergens.github.io/RdistNN/articles/conditional-z-estimation.md);
return here when cross-estimator names or inferential scales become
ambiguous.
