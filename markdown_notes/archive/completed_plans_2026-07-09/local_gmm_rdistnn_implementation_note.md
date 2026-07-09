> **Status (2026-07-09): implemented and archived.** `gmm_dnn()` and
> `gmm_model_dr_iv()` exist as experimental companion APIs, following the
> estimator-layer design in this note.

# Local GMM for `RdistNN`: implementation note

Date: 2026-03-30

## Bottom line

The clean way to add the local GMM variant from `Cond_GMM_with_DNN` to `RdistNN` is to add a new estimator layer above the existing DNN weight engine, not to treat it as a small extension of `plqr_dnn()` or of the current scalar `score` API.

`RdistNN` already has most of the localization machinery we need:

- `weights_matrix()` / `weights_dnn()` for DNN weights
- batched target handling through `X0`
- reusable ranks / orders / distances
- `n_eff = 1 / sum(W^2)` style diagnostics
- delete-1 rank-update machinery that could later support a full refit jackknife

What it does not yet have is the GMM-specific layer:

- a vector-moment model abstraction
- a two-step optimizer
- local Jacobian and local covariance estimators
- plug-in sandwich inference
- overidentification diagnostics such as the local `J`-statistic

My preferred first version is:

- keep the current normalized-weight convention used everywhere else in the package
- accept nuisance inputs from the user rather than fitting them internally
- implement the paper's identity first step plus efficient second step
- return sandwich SEs and the local `J`-test
- treat delete-1 jackknife for the final estimator as a second-phase feature, not the phase-1 default

## Relevant current files

The package files that matter most for this extension are:

- `R/api-weights.R`
- `R/user-dnn-core.R`
- `R/user-dnn-plqr.R`
- `R/score-plqr.R`
- `vignettes/extending.Rmd`
- `vignettes/conditional-z-estimation.Rmd`

On the theory side, the main source is:

- `../Cond_GMM_with_DNN/Paper/Parts/02_Setup.tex`
- `../Cond_GMM_with_DNN/Paper/Parts/03_Results.tex`

## What the current package architecture already gives us

### 1. The weight layer is already in the right place

`weights_matrix()` and `weights_dnn()` already separate:

- ranking / distance evaluation
- rank-to-weight mapping
- multi-target batching

That is exactly the right abstraction for local GMM. Nothing about the GMM paper suggests rewriting the DNN weight layer.

### 2. The current user wrappers are fixed-signal linear aggregators

The current high-level wrappers all have the same computational shape:

- build an observation-level signal matrix
- compute `W`
- return `crossprod(signal, W)`

This is true for:

- `np_dnn()`
- `cate_dnn()`
- `plqr_dnn()`

The fast SE helper in `R/user-dnn-core.R`, `.dnn_jackknife_se_loo_matrix()`, is built for exactly that pattern: a fixed signal matrix is re-aggregated under delete-1 reweighting.

That matters because local GMM is not of that form once the parameter is estimated by solving an optimization problem.

### 3. The score layer is useful, but not the right top-level abstraction

The score API is modular and well designed, but it is still scalar-output and subset-based.

The built-in pure-score factories such as `score_cate_rho()` and `score_plqr_component()`:

- stash per-observation contributions `g_i` in `setup()`
- return subset sums in `eval()` / `eval_batch()`
- leave all weighting outside the score object

That is useful, but local GMM needs a higher-level object that naturally returns:

- an `n x q` observation-level moment matrix
- optionally an observation-level Jacobian contribution
- model-specific starting values / constraints

Trying to force `q > 1` local GMM into "one scalar score object per moment component" would be possible, but it would be awkward for:

- Jacobian handling
- `Omega_hat`
- optimizer code
- local `J`-statistics
- multi-moment diagnostics

My conclusion is that the score layer is a helper here, not the main abstraction.

## What the GMM repo implies for code

The GMM paper defines the local estimator through the weighted conditional moment average

- `gbar(theta; x0)`

followed by:

1. an identity-weighted first step
2. a local covariance estimate `Omega_hat`
3. an efficient second step using `Omega_hat^{-1}`

Translated into package code, the core quantities should be:

- `g_obs(theta)`: `n x q` observation-level moment matrix
- `gbar(theta) = crossprod(g_obs(theta), W)` under the package's normalized `W`
- `Omega_hat(theta) = sum_i W_i g_i(theta) g_i(theta)'`
- `D_hat(theta) = d gbar(theta) / d theta`
- `V_hat = (D_hat' Omega_hat^{-1} D_hat)^{-1}`
- `J = n_eff * gbar(theta_hat)' Omega_hat^{-1} gbar(theta_hat)`

That is the right computational translation of the paper for `RdistNN`.

## A notation issue that should be resolved explicitly in code

The GMM paper setup currently writes some quantities with an extra `1 / n` factor while also describing the DNN weights as living in the simplex and using

- `n_eff = (sum_i w_i^2)^(-1)`

For package implementation I would not try to mirror that notation literally. I would follow the existing `RdistNN` convention:

- columns of `W` sum to one
- all local averages are computed as weighted sums `crossprod(., W)`

If we do not make that translation explicit, it is easy to scale `gbar`, `Omega_hat`, `V_hat`, and the local `J`-statistic incorrectly.

## Main implementation problems I see

### 1. Local GMM is nonlinear; the current fast jackknife helper is not enough

`plqr_dnn()` computes moments at a user-supplied `beta`. It does not solve for `beta`.

By contrast, the local GMM estimator changes because:

- the criterion is optimized
- the second-step weighting matrix depends on the first-step estimate
- the final estimator is not a linear map of a fixed signal matrix

So I would not reuse `.dnn_jackknife_se_loo_matrix()` as the phase-1 inference path. That helper is right for localized averages of a fixed signal, but not for a refit two-step GMM estimator.

The correct phase-1 default is the plug-in sandwich estimator from the paper.

If we later want jackknife inference, it should be a true refit jackknife of the whole estimator.

### 2. We need a vector-moment model contract

The right object is not "a scalar score," but something like:

- `moment_obs(data, theta, ctx) -> n x q`
- `jacobian_obs(data, theta, ctx) -> n x (q x p)` or array

plus optional:

- `setup(data, ...)`
- `start(data, ctx, x0)`
- parameter bounds / projection
- names for coefficients and moments

Without that contract, the optimizer and the inference code will end up re-implementing model-specific logic inside the wrapper.

### 3. The nuisance contract has to be explicit

The paper is written with leave-one-out nuisance evaluations `eta_hat_{-i}` inside the local moments.

`RdistNN` currently does not manage LOO nuisance fitting internally. Instead, the wrappers mostly consume nuisance values that are already available as columns or vectors.

I think local GMM should follow the same phase-1 pattern:

- the user supplies nuisance evaluations already aligned with observations
- the package treats them as data inputs

That is especially important because different GMM models will want very different nuisance learners.

If we later add package-managed nuisance fitting, that should be a separate layer.

### 4. Jacobian estimation needs an analytic-first / FD-fallback path

Best case:

- the model supplies observation-level Jacobian contributions, so `D_hat` is just a weighted average

Fallback:

- finite-difference the local moment map `gbar(theta)`, not each individual `g_i(theta)`

That second option is important. It lets us support generic models without requiring an `n x q x p` Jacobian array.

So the implementation contract should be:

- analytic Jacobian if available
- finite-difference `gbar(theta)` otherwise

This is simpler than forcing every model to provide observation-level derivatives.

### 5. `Omega_hat` can be locally singular or nearly singular

The overidentified setting creates numerical issues that the current wrappers do not have.

We will need:

- explicit symmetrization of `Omega_hat`
- positive-definiteness checks
- eigenvalue flooring or ridge regularization
- clear warnings when the local moment covariance is nearly singular

The paper's many-moments condition

- `q^2 / n_eff -> 0`

should also become a runtime diagnostic. Even if the asymptotics are fine on paper, local `q / n_eff` can be practically large.

### 6. Multi-target support will need warm starts

`weights_dnn()` already supports many targets at once. A local GMM loop over many `x0` values can reuse that, but solving two optimizations per target will still be expensive.

The obvious phase-1 approach is:

- compute weights for all targets
- loop over targets
- use the previous target's solution as the next starting value

That should be especially effective when `X0` is a grid.

### 7. The result object needs to be richer than the existing wrappers

A local GMM result should contain at least:

- first-step estimate
- second-step estimate
- `Omega_hat`
- `D_hat`
- `V_hat`
- `se`
- local `J`
- `df = q - p`
- `n_eff`
- optimizer convergence info
- optional `W`, `ranks`, `ord`, `dx`

This is a larger object than `np_result`, `cate_result`, or `plqr_result`, but that is appropriate because the estimator is more complex.

## Why I would stay on the L-stat side computationally

The current GMM paper is theoretically connected to the generalized-U-statistic machinery, but the actual computation should stay on the observation-level weighted-average side.

That is already the package's preferred implementation style:

- compute observation-level contributions
- aggregate them with the L-stat weight matrix `W`

I would not try to compute local GMM by explicit subset enumeration or by building a high-level solver around subset kernels. The score factories are useful for diagnostics and experimentation, but the production implementation should stay on the same L-stat side as the existing user wrappers.

## A concrete phase-1 target: local DR IV with multiple instruments

The DR-IV example in the GMM paper is the best first implementation target.

Suppose the user supplies:

- outcome `Y`
- endogenous treatment `D`
- instrument columns `Z_1, ..., Z_q`
- nuisance columns `pi_1, ..., pi_q`
- nuisance column `muY`
- nuisance column `muD`

For scalar `theta`, the observation-level moment vector is

- `g_i(theta) = (Z_i - pi_i) * (Y_i - muY_i - (D_i - muD_i) * theta)`

where the multiplication is componentwise over the `q` instruments.

The Jacobian is then

- `G_i(theta) = -(Z_i - pi_i) * (D_i - muD_i)`

which is especially convenient because it does not depend on `theta`.

This gives us:

- a clean overidentified example
- an immediate local `J`-statistic
- a natural `q > 1` benchmark
- an easy reduction to the just-identified `q = 1` case

It is a much better first smoke test than trying to start with a fully generic black-box GMM model.

## Design I would recommend

### Core abstraction: a model object

I would add something like `new_gmm_model()` with fields:

- `name`
- `setup(data, ...)`
- `moment_obs(data, theta, ctx)`
- `jacobian_obs(data, theta, ctx)` optional
- `start(data, ctx, x0)` optional
- `project(theta)` or bounds optional
- `coef_names`
- `moment_names`

This keeps the package layered cleanly:

- weights stay in the weight layer
- model-specific moments stay in the model layer
- optimization and inference stay in the estimator layer

### Core estimator: `gmm_dnn()` or `local_gmm_dnn()`

Arguments I would expose first:

- `data`, `x_cols`, `X`, `X0`, `x0`
- `model`
- `distance`, `method`, `k`, `s`, `normalize`
- `start`
- `optimizer = c("BFGS", "Nelder-Mead", "L-BFGS-B")`
- `jacobian = c("analytic", "fd")`
- `omega_regularization`
- `se = c("sandwich", "none")`
- `return = c("W", "ranks", "ord", "dx", "g_obs", "omega", "jacobian", "ess")`

The engine should:

1. build `W` once per target
2. solve the identity-weighted first step
3. form `Omega_hat` at the first-step solution
4. solve the efficient second step
5. build `D_hat`, `V_hat`, `se`, and `J`

### Inference default

I would not make jackknife the default for local GMM.

The current paper's feasible inference is sandwich-based, and the current fast jackknife path is designed for fixed observation-level signals. A full delete-1 refit jackknife could be valuable later, but it is a distinct feature.

### Weight convention

I would keep the current package convention:

- `colSums(W) = 1`
- `n_eff = 1 / colSums(W^2)`
- all local averages are weighted sums using `W`

That makes local GMM consistent with the rest of `RdistNN`.

## Testing strategy I would want

I would add tests in this order:

1. A toy linear-GMM example where the optimizer has a known solution.
2. A just-identified case `q = p` where first-step and second-step estimates agree.
3. A case where `Omega_hat` is proportional to the identity, so the second step collapses to the first.
4. Analytic-vs-FD Jacobian parity tests.
5. Invariance under permutation of the moment components.
6. Regularization tests for nearly singular `Omega_hat`.
7. A local IV example with `q = 1`, `p = 1` that matches the weighted Wald ratio.
8. Multi-target tests verifying that warm starts do not change the final answer.

## Recommended implementation order

1. Settle the API and object names before writing code.
2. Implement a model object plus a single-target `gmm_dnn()` prototype.
3. Make the first version accept user-supplied nuisance inputs only.
4. Add analytic Jacobian support plus FD fallback on `gbar(theta)`.
5. Add sandwich SEs, `vcov()`, `confint()`, and the local `J`-statistic.
6. Add multi-target looping with warm starts.
7. Add a DR-IV factory / demo dataset.
8. Only after that, consider refit jackknife, CUE, or built-in nuisance learners.

## My preferred path

If I had to choose one concrete implementation path right now, it would be:

- leave the current DNN weight engine untouched
- add a new `gmm_dnn()` layer instead of mutating `plqr_dnn()`
- use a dedicated vector-moment model object instead of forcing everything through scalar score objects
- make sandwich inference the default
- start with externally supplied nuisance inputs and the paper's DR-IV example
- leave refit jackknife and automatic nuisance fitting for a second pass

That path is closest to the GMM paper, fits the existing `RdistNN` architecture, and keeps the first implementation small enough to land cleanly.

## Implementation decisions for the first estimator layer

The first implementation layer uses the following public names and defaults:

- `new_gmm_model()` / `is_gmm_model()` define the vector-moment extension point.
- `gmm_model_dr_iv()` is the first concrete model factory.
- `gmm_dnn()` is the user-facing estimator wrapper.
- Estimator outputs use generic `theta` naming. IV simulation DGPs may keep
  `tau_x0` as an alias, but simulation metrics and estimator results use
  `theta`, `theta_hat`, and `theta0`.
- `gmm_dnn()` consumes nuisance values supplied in `data` or `model_data`; it
  does not fit nuisances internally.
- Sandwich inference is the first supported inference path. Delete-1 refit
  jackknife, CUE, and automatic nuisance learning are deferred.
- Local moment covariance matrices are symmetrized and ridge-regularized with
  `omega_ridge = 1e-8` by default before inversion.
