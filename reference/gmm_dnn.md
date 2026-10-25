# Local GMM estimation via DNN rank weights

Estimates a local parameter by solving a two-step GMM problem after
localizing observations with DNN rank weights. The first step uses an
identity moment weighting matrix; the second step uses the local moment
covariance estimated at the first-step solution. This local IV-GMM
interface is experimental and may change.

## Usage

``` r
gmm_dnn(
  data = NULL,
  x_cols = NULL,
  X = NULL,
  X0 = NULL,
  x0 = NULL,
  model,
  model_data = NULL,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L,
  s = NULL,
  tuning = NULL,
  normalize = TRUE,
  start = NULL,
  optimizer = c("BFGS", "Nelder-Mead", "L-BFGS-B"),
  jacobian = c("auto", "analytic", "fd"),
  omega_ridge = 1e-08,
  fd_step = 1e-05,
  se = TRUE,
  ctx = NULL,
  return = c("W")
)
```

## Arguments

- data:

  Optional data frame containing covariates and model variables. Column
  names must be non-missing, non-empty, and unique.

- x_cols:

  Names or indices of covariate columns in `data`.

- X:

  Numeric `n x p` covariate matrix for the matrix interface.

- X0:

  Numeric matrix of target points; mutually exclusive with `x0`.

- x0:

  Numeric vector for one target point.

- model:

  A
  [`new_gmm_model()`](https://JakobJuergens.github.io/RdistNN/reference/new_gmm_model.md)
  object.

- model_data:

  Optional model data used with the matrix interface.

- distance:

  Distance name or specification.

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
  omitted.

- normalize:

  Logical; normalize weight columns to sum to one.

- start:

  Optional numeric vector of length `p_theta`, or `p_theta x m` matrix,
  used as optimizer starts. When omitted, the model's
  [`start()`](https://rdrr.io/r/stats/start.html) callback is evaluated
  separately at every target, or zero is used when the model has no
  callback. Targets do not inherit previous targets' starts or
  solutions. An explicit vector is reused at all targets; matrix column
  `j` supplies target `j`. Either explicit form takes precedence over
  the model callback.

- optimizer:

  Optimizer method passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html). Both stages
  allow 1000 iterations. BFGS and Nelder-Mead use `reltol = 1e-12`;
  L-BFGS-B uses `factr = 1e-12 / .Machine$double.eps` and the default
  projected-gradient tolerance. Convergence statuses remain available in
  the returned object.

- jacobian:

  One of `"auto"`, `"analytic"`, or `"fd"`.

- omega_ridge:

  Non-negative ridge multiplier for local moment covariance.

- fd_step:

  Positive finite finite-difference step multiplier.

- se:

  Logical; compute sandwich standard errors.

- ctx:

  Optional context passed to distance and ranking components.

- return:

  Optional extras:
  `c("W", "ranks", "ord", "dx", "g_obs", "omega", "jacobian", "ess", "neighbors", "tuning")`.
  `neighbors` is a `k x m` matrix of nearest training-observation
  indices, not the full DNN weight support.

## Value

A `rdistnn_gmm_result` object with estimates, diagnostics, metadata, the
`m x p` evaluation-point matrix `x0` (with covariate column names
propagated from the training data when available), and requested extras.

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
Refit saved models after package updates during experimental
development.

Two interfaces are supported.

1.  Data-frame interface: pass `data`, `x_cols`, and a
    `rdistnn_gmm_model`.

2.  Matrix interface: pass `X`, one of `X0` or `x0`, a
    `rdistnn_gmm_model`, and `model_data` when the model needs
    observation-level variables.

Nuisance quantities are not fit internally. They must be supplied to the
model through `data` or `model_data`.

The resolved distance specification and its fitted training context are
retained for weight reconstruction. Replacing a registered metric after
fitting does not change the fit's geometry. Custom callbacks with
mutable external state retain their own reproducibility obligations.

Model callbacks must return their documented finite numeric values at
the supplied starts and optimizer trial points. Callback errors,
malformed returned values, interrupts and caller-promoted warnings stop
the fit. Only nonfinite quadratic arithmetic from otherwise valid
moments receives an optimizer penalty. A local moment covariance that
cannot be inverted leaves that target unavailable, with its reason in
`convergence` and one partial-result warning per call. Singular local
information retains the coefficient estimate with unavailable
uncertainty; information outside the numeric range also gives a
rescaling warning. No callback-domain recovery interface or automatic
retry at a different start is supplied.

Requested Jacobian diagnostics can retain arithmetic nonfinite values
from differencing finite callback evaluations when coefficient
uncertainty is unavailable or deliberately disabled. This does not relax
the finite-value contract for callback returns.

Exact DNN localization uses the full rank permutation at each target
rather than an `s`-nearest-neighbor truncation. For exact DNN fits, the
unaveraged rank-table support edge is rank `n - s + k`, not rank `s`.

Averaging can extend nonzero observation weights through the end of a
distance block that straddles the nominal rank-table edge. Effective
sample sizes and weight diagnostics use the resulting observation
weights.

## Examples

``` r
dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 120, kx = 2, q = 2, seed = 1)
model <- gmm_model_dr_iv(
  z_cols = c("Z1", "Z2"),
  pi_cols = c("piZ1_oracle", "piZ2_oracle"),
  muY_col = "muY_oracle",
  muD_col = "muD_oracle"
)
fit <- gmm_dnn(
  data = dgp$data,
  x_cols = c("X1", "X2"),
  x0 = dgp$x0,
  model = model,
  k = 1, s = 10,
  return = "ess"
)
fit$theta
#>       target_1
#> theta 1.196561
```
