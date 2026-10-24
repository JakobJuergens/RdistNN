# CATE via DNN rank-weighted pseudo-outcome averaging

Estimate conditional average treatment effects by combining a
pseudo-outcome with DNN-based local weights.

## Usage

``` r
cate_dnn(
  data = NULL,
  x_cols = NULL,
  y_col = NULL,
  t_col = NULL,
  m0_col = NULL,
  m1_col = NULL,
  e_col = NULL,
  pseudo_outcome = NULL,
  X = NULL,
  X0 = NULL,
  x0 = NULL,
  clip = 1e-08,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L,
  s = NULL,
  tuning = NULL,
  tuning_mode = c("default", "grid"),
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- data:

  Data frame with covariates and optional nuisance columns. Column names
  must be non-missing, non-empty, and unique.

- x_cols:

  Names or indices of covariate columns in `data`.

- y_col:

  Name or index of the outcome column in `data`.

- t_col:

  Name or index of the treatment column in `data` (must be 0/1).

- m0_col:

  Name or index of the supplied control-arm conditional outcome mean.

- m1_col:

  Name or index of the supplied treated-arm conditional outcome mean.

- e_col:

  Name or index of the supplied treatment propensity predictions.

- pseudo_outcome:

  Numeric vector of length `n`. If supplied, nuisance columns are not
  required and are ignored. The signal is localized unchanged; its
  construction and inferential justification are the caller's
  responsibility.

- X:

  Numeric matrix (`n x p`) for the matrix interface.

- X0:

  Numeric matrix (`m x p`) of target points; mutually exclusive with
  `x0`.

- x0:

  Numeric vector of length `p` for a single target point.

- clip:

  Finite propensity clipping constant in `(0, 0.5)` with `1 - clip < 1`
  in double precision, default `1e-8`. Raw nuisance propensities must be
  finite and lie in `[0, 1]`; nuisance-based AIPW scores use
  `pmin(pmax(e, clip), 1 - clip)`. This argument does not alter a
  supplied `pseudo_outcome`.

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
  only. `"grid"` evaluates the feasible candidates from `tuning`, or
  from
  [`tuning_dnn_paper()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_dnn_paper.md)
  when `tuning` is omitted, and stores candidate-level estimates,
  standard errors, diagnostics, and metadata in `tuning_grid`. Grid mode
  cannot be combined with explicit `s`.

- normalize:

  Logical; renormalize columns of `W` to sum to 1.

- se:

  Logical; compute delete-1 jackknife standard errors.

- jackknife:

  One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default) uses cached
  full-sample ranks with delete-1 rank updates, `"loo"` recomputes
  leave-one-out weights directly, and `"none"` skips SEs.

- ctx:

  Optional context passed to distance and ranking components.

- return:

  Optional extras in addition to the always-returned `tau`, `se`, and
  `meta`. Allowed values:
  `c("W", "ranks", "ord", "dx", "psi", "nuisance", "ess", "neighbors", "tuning")`.
  `neighbors` is a `k x m` matrix of nearest training-observation
  indices, one column per target. It is a rank summary, not the full DNN
  weight support.

## Value

An object of class `rdistnn_cate_result` with components:

- `tau`: estimated CATE values at target points;

- `se`: optional jackknife standard errors;

- `x0`: the `m x p` matrix of evaluation points, with covariate column
  names propagated from the training data when available;

- `meta`: metadata (`n`, `m`, `p`, `k`, `s`, `method`, `distance`,
  `signal_source`, `created_with_version`); see Signal provenance;

- `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
  `diagnostics`, and `meta` tables. Candidate diagnostics report top-5
  weight concentration separately from the estimator's `k`;

- `.model`: retained observations, supplied signals and fitted geometry
  used by prediction and joint covariance, independently of optional
  extras;

- optional extras requested via `return`.

## Details

The documented sampling framework uses independent observations; cluster
and survey-design inference are not supplied. Observations with exactly
equal computed distance to a target share their block's average rank
weight, even when their covariates differ. Groups are target-specific
and use exact represented equality, without a tolerance. Ordering uses
distance then input row index; outcomes and nuisances do not define
groups. Every jackknife deletion averages the surviving blocks with the
sample-size-n-1 rule. A warning is emitted once per public call with
distance ties: the accompanying theory assumes a continuous distance
distribution and does not establish inference under persistent ties.
Earlier retained fit schemas require refitting.

Two interfaces are supported.

1.  Data-frame interface (recommended): pass `data` and `x_cols`, then
    provide either `pseudo_outcome` directly or nuisance columns used to
    construct it. For AIPW nuisance construction, use `y_col`, `t_col`,
    `e_col`, `m0_col`, and `m1_col`.

2.  Matrix interface: pass `X`, one of `X0` or `x0`, and
    `pseudo_outcome`.

As in
[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
exact DNN localization uses the full rank permutation at each target
rather than an `s`-nearest-neighbor truncation. The `"fast"` jackknife
option reuses the full-sample rank matrix and updates ranks under
deletion. For exact DNN fits, the unaveraged rank-table support edge is
rank `n - s + k`, not rank `s`.

Averaging can extend nonzero observation weights through the end of a
distance block that straddles the nominal rank-table edge. Effective
sample sizes and weight diagnostics use the resulting observation
weights.

`tuning_mode = "default"` evaluates only the default kernel order and
preserves the ordinary scalar-fit return shape. `tuning_mode = "grid"`
evaluates the feasible candidates in `tuning`, or the paper default
tuning grid when `tuning` is omitted, and attaches a tidy `tuning_grid`
payload. The pseudo-outcome and nuisance columns are resolved once
before grid evaluation. The primary `tau`, `se`, `meta`, and optional
extras such as `W`, `ranks`, `ess`, `psi`, and `nuisance` remain
anchored to the default candidate. Full per-candidate weight matrices
are not stored. Candidate diagnostics include `ess`, `max_weight`,
`top_n`, `top5_weight_share`, and `degenerate`; `top_n` is fixed at 5
and is independent of the estimator's `k`.

For causal identification and pointwise inference, including conditions
on nuisance accuracy and deletion stability, see the accompanying
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
This function accepts supplied oracle, out-of-fold or full-sample
predictions; it does not verify those conditions for a learner. The
normalized Robinson construction is not implemented; use AIPW inputs or
a justified supplied signal. A pooled outcome mean cannot replace the
two arm-specific means.

With `data`, omit `X`; a direct `pseudo_outcome` vector remains
supported. Column selectors require `data`. Prediction retains the
fitted metric and context even if a registered distance name is
subsequently replaced.

## Signal provenance

`meta$signal_source` records the actual input branch: `"aipw"` for
built-in nuisance construction, `"supplied"` for `pseudo_outcome` (which
takes precedence over nuisance selectors). Neither signal values nor
column names establish provenance. `meta$created_with_version` records
the package version loaded when the original fit was created. Summaries
and candidate views preserve these fields. Provenance does not certify
causal identification, nuisance accuracy/stability, localization bias,
or pointwise inference.

## See also

[`weights_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/weights_dnn.md),
[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md),
[`score_cate_rho()`](https://JakobJuergens.github.io/RdistNN/reference/score_cate_rho.md),
[`cate_demo_data()`](https://JakobJuergens.github.io/RdistNN/reference/cate_demo_data.md),
[`vcov.rdistnn_cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/vcov.rdistnn_cate_result.md),
[`predict.rdistnn_cate_result()`](https://JakobJuergens.github.io/RdistNN/reference/predict.rdistnn_cate_result.md).

## Examples

``` r
set.seed(1)
n <- 50
x_mat <- matrix(rnorm(n * 2), n, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
psi <- rnorm(n)
out <- cate_dnn(
  X = x_mat, X0 = x0_mat, pseudo_outcome = psi,
  distance = "l2", method = "exact", k = 1L, s = 3L,
  se = FALSE, return = c("W", "ess")
)
out$tau
#> [1] -0.02293282 -0.04260484 -0.24413829
```
