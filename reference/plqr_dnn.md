# PLQR local moment evaluation and direct coefficient inference via DNN rank weights (experimental)

Compute localized orthogonal PLQR moments \$\$ \mathbb{E}\\\left\[
\left(\tau - \mathbf{1}\\Y \le D^\top \theta + g_0(U)\\\right) \left(D -
m_D(U)\right) \\\middle\|\\ X \approx x_0 \right\] \$\$ at a supplied
coefficient value `theta_eval`, using DNN rank weights at one or
multiple target points. When positive threshold density values are
supplied through `density` or `density_col`, the same fit reports the
direct PLQR coefficient sandwich used for standard errors and confidence
intervals. No local root solving is performed: `theta` in the result is
the supplied evaluation value (`meta$theta_source = "supplied"`).

## Usage

``` r
plqr_dnn(
  data = NULL,
  x_cols = NULL,
  D_cols = NULL,
  g_col = NULL,
  mD_cols = NULL,
  y_col = NULL,
  theta_eval,
  tau = 0.5,
  X = NULL,
  D = NULL,
  g = NULL,
  mD = NULL,
  Y = NULL,
  density = NULL,
  density_col = NULL,
  X0 = NULL,
  x0 = NULL,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L,
  s = NULL,
  tuning = NULL,
  tuning_mode = c("default", "grid"),
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  theta_meat = c("analytic", "jackknife"),
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- data:

  Optional data frame containing covariates, regressors, outcome, and
  nuisance columns. Column names must be non-missing, non-empty, and
  unique.

- x_cols:

  Names or indices of covariate columns in `data`.

- D_cols:

  Names or indices of partially linear regressor columns in `data`.

- g_col:

  Name or index of the shape-nuisance column in `data` holding the
  values of \\g_0(U)\\.

- mD_cols:

  Names or indices of nuisance columns for `m_D(U)`.

- y_col:

  Name or index of the outcome column in `data`.

- theta_eval:

  Numeric coefficient vector at which the PLQR moment is evaluated, with
  length `length(D_cols)` (or `ncol(D)` in the matrix interface). No
  root solving is performed; the result's `theta` repeats this supplied
  value across targets.

- tau:

  Quantile level in `(0,1)`.

- X:

  Numeric matrix (`n x p`) of covariates for the matrix interface.

- D:

  Numeric matrix (`n x p_d`) of partially linear regressors for the
  matrix interface.

- g:

  Numeric vector (`n`) with \\g_0(U)\\ values for the matrix interface.

- mD:

  Numeric matrix (`n x p_d`) with `m_D(U)` values for the matrix
  interface.

- Y:

  Numeric outcome vector of length `n` for the matrix interface.

- density:

  Optional positive numeric vector (`n`) estimating the conditional
  density of `Y` at the threshold `D' theta_eval + g`. When supplied
  with `se = TRUE`, PLQR uses the direct coefficient sandwich route and
  returns coefficient standard errors in `theta_se`.

- density_col:

  Optional column name or index in `data` containing the positive
  density values used for the direct coefficient sandwich.

- X0:

  Numeric matrix (`m x p`) of target points; mutually exclusive with
  `x0`.

- x0:

  Numeric vector of length `p` for a single target point.

- distance:

  Distance name or specification (for example `"l2"`).

- method:

  One of `"exact"` or `"asymp"` for DNN weights.

- k:

  Integer DNN neighbour count inside each size-`s` subsample. The paper
  theory uses `k = 1` (the default).

- s:

  Integer kernel order (undersmoothing parameter), or `NULL`. When both
  `s` and `tuning` are `NULL` (the default), `s` is resolved from the
  paper undersmoothing rule
  [`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md),
  \\s = \lfloor 0.4 n^{(k_x+1)/(k_x+2)} \rfloor\\. Note that small fixed
  values such as `s = 1` (or `k = s`) yield uniform weights and hence no
  localization; a warning is emitted in that case.

- tuning:

  Optional kernel-order tuning specification. When supplied, `s` must be
  omitted and is resolved from `tuning`.

- tuning_mode:

  One of `"default"` or `"grid"`. `"default"` evaluates the default `s`
  only. `"grid"` evaluates the feasible candidates from `tuning` (or the
  paper default tuning grid when `tuning` is `NULL`) and stores
  candidate-level moments, standard errors, diagnostics, and metadata in
  `tuning_grid`. Grid mode cannot be combined with explicit `s`.

- normalize:

  Logical; renormalize columns of `W` to sum to 1.

- se:

  Logical; compute delete-1 jackknife standard errors.

- jackknife:

  One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default) uses cached
  full-sample ranks with delete-1 rank updates, `"loo"` recomputes
  leave-one-out weights directly, and `"none"` skips SEs.

- theta_meat:

  Meat estimator for the coefficient sandwich when
  `density`/`density_col` is supplied. `"analytic"` (default) uses the
  analytic score covariance \\S(x) = \tau(1-\tau)\sum_i w_i \tilde D_i
  \tilde D_i'\\; `"jackknife"` uses the full delete-1 jackknife
  covariance of the localized vector moment (`moment_vcov`) as the meat,
  so \\\widehat V\_\theta = J^{-1} \widehat S\_{JK} J^{-T}\\. The bread
  is density-backed in both routes. `"jackknife"` requires
  `density`/`density_col`, `se = TRUE`, and `jackknife != "none"`.

- ctx:

  Optional context passed to distance and ranking components.

- return:

  Optional extras in addition to the always-returned `moment`, `theta`,
  `theta_eval`, `tau`, and `meta`. Allowed values:
  `c("W", "ranks", "ord", "dx", "score", "ess", "neighbors", "tuning")`.
  The `neighbors` extra returns the first `k` nearest observations used
  inside each size-`s` subsample, not the full DNN weight support.

## Value

An object of class `plqr_result` with components:

- `moment`: matrix (`p_d x m`) of localized PLQR moments;

- `theta`: matrix (`p_d x m`) of the supplied coefficient values,
  repeated across targets until a root-solving path is requested in
  future versions;

- `theta_eval`: the supplied coefficient vector used in the moment map;

- `theta_se`: coefficient standard errors (`p_d x m`), present only when
  `density` or `density_col` supplies the direct PLQR sandwich route;

- `moment_se`: componentwise jackknife standard errors for `moment`.
  Present when no density is supplied and `se = TRUE` with
  `jackknife != "none"` (moment-level diagnostic only), and when
  `theta_meat = "jackknife"` (as the diagonal of `moment_vcov`). A
  componentwise `moment_se` is never by itself coefficient inference;

- `moment_vcov`: with `theta_meat = "jackknife"`, the full delete-1
  jackknife covariance of the localized vector moment, one `p_d x p_d`
  matrix per target — the sandwich meat;

- `tau`: quantile level;

- `x0`: the `m x p` matrix of evaluation points, with covariate column
  names propagated from the training data when available;

- `meta`: metadata (`n`, `m`, `p`, `p_d`, `k`, `s`, `method`,
  `distance`, `theta_source = "supplied"`, `theta_se_method`,
  `theta_meat_method`, `moment_se_method`);

- when direct coefficient inference is available, `theta_vcov`,
  `theta_omega`, `jacobian`, and `neff` (the paper's effective sample
  size \\n\_{\mathrm{eff}}(x) = 1/\sum_i w_i(x)^2\\, identical to the
  `ess` return extra), plus `score_covariance` (the analytic meat) when
  `theta_meat = "analytic"`;

- `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
  `diagnostics`, and `meta` tables. Candidate diagnostics report top-5
  weight concentration separately from the estimator's `k`;

- optional extras requested via `return`.

## Details

The documented sampling framework uses independent observations; cluster
and survey-design inference are not supplied. Input row order breaks
distance ties, so permuting tied observations can change estimates.

Two interfaces are supported.

1.  Data-frame interface (recommended): pass `data` with `x_cols`,
    `D_cols`, `g_col`, `mD_cols`, and `y_col`.

2.  Matrix interface: pass `X`, `D`, `g`, `mD`, and `Y`.

No nuisance models are fit internally; this wrapper assumes the
partially-linear shape nuisance \\g_0(U)\\ and the residualizing
nuisance \\m_D(U) = E\[D \mid U\]\\ are supplied. The threshold is
\\D'\theta + g_0(U)\\; \\g_0\\ alone is not itself a quantile function.
Both `g` and `mD` may depend on covariates beyond `x_cols` (the paper's
\\U = (X, Z)\\); the package does not restrict their domain, since they
are supplied pre-computed. Density values for direct coefficient
inference are also supplied by the user rather than estimated
internally; they play the role of the paper's local conditional density
at the threshold, entering only the coefficient Jacobian. See
`markdown_notes/design/15_plqr_variance_estimation_rationale.md` for the
full nuisance-notation mapping and for why coefficient-level inference
requires this direct-sandwich route rather than the jackknife SE used by
[`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)/[`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md).
Direct coefficient inference requires `normalize = TRUE`, matching the
theorem convention of weights summing to one per target.

With a density input, the coefficient sandwich supports two meat
estimators (`theta_meat`): the analytic score covariance (default) and
the full delete-1 jackknife covariance of the localized vector moment.
The bread is density-backed in both routes; only the meat changes. A
componentwise `moment_se` alone is never coefficient inference — see
`markdown_notes/design/15_plqr_variance_estimation_rationale.md` for the
route table.

As in
[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
exact DNN localization uses the full rank permutation at each target
rather than an `s`-nearest-neighbor truncation. The `"fast"` jackknife
option reuses the full-sample rank matrix and updates ranks under
deletion. For exact DNN fits, the finite-population support edge is rank
`n - s + k`, not rank `s`.

`tuning_mode = "default"` evaluates only the default kernel order and
preserves the ordinary scalar-fit return shape. `tuning_mode = "grid"`
evaluates the feasible candidates in `tuning`, or the paper default
tuning grid when `tuning` is omitted, and attaches a tidy `tuning_grid`
payload. The payload stores one long-form row per candidate, target, and
PLQR moment component; the public `moment` and `moment_se` matrix fields
remain anchored to the default candidate. Candidate diagnostics include
`ess`, `max_weight`, `top_n`, `top5_weight_share`, and `degenerate`;
`top_n` is fixed at 5 and is independent of the estimator's `k`. Direct
PLQR sandwich inference from `density` or `density_col` currently uses
`tuning_mode = "default"` only.

## Lifecycle

**Experimental.** PLQR functionality is outside the CATE and regression
workflow focus of the initial public release. Its API and inference
interfaces may change in later releases.

## See also

[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
[`plqr_score_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_score_matrix.md),
[`plqr_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_demo_data.md).

## Examples

``` r
set.seed(1)
demo <- plqr_demo_data(n = 120, p = 2, p_d = 2, tau = 0.5, seed = 1)
out <- plqr_dnn(
  data = demo$data,
  x_cols = c("X1", "X2"),
  D_cols = c("D1", "D2"),
  g_col = "g_oracle",
  mD_cols = c("mD1", "mD2"),
  y_col = "Y",
  density_col = "density_oracle",
  theta_eval = demo$truth$theta_true,
  tau = 0.5,
  x0 = c(0, 0),
  method = "exact",
  k = 1L, s = 3L,
  se = TRUE
)
summary(out)
#> <summary.plqr_result>
#>   targets: 1
#>   theta_source: supplied; no local PLQR root solving was performed
#>     target term estimate std.error statistic p.value conf.low conf.high
#> 1 target_1   D1      0.5   0.04769        NA      NA   0.4065    0.5935
#> 2 target_1   D2      0.5   0.05041        NA      NA   0.4012    0.5988
```
