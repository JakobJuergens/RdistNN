# Run the CATE simulation grid

Expands a grid of DGP designs, sample sizes, and covariate dimensions,
then runs `n_reps` independent Monte Carlo replications for each cell
using
[`sim_cate_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_one_rep.md).
Each replicate draws one dataset and evaluates **all** requested
`methods` on the same draw, enabling like-for-like comparisons. Returns
a tidy data.frame of raw per-replicate, per-method results. Use
[`summarise_sim_results()`](https://JakobJuergens.github.io/RdistNN/reference/summarise_sim_results.md)
to aggregate to one row per
`(method, design, n, kx, candidate_id, lambda, s)`.

## Usage

``` r
sim_cate_run(
  designs = c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap", "dgp4_clustered",
    "dgp5_correlated"),
  n_vals = c(500L, 1000L, 2000L),
  kx_vals = 2L,
  n_reps = 500L,
  master_seed = 42L,
  k = 1L,
  s_fn = function(n, kx) floor(0.4 * n^((kx + 1L)/(kx + 2L))),
  tuning = NULL,
  tuning_mode = c("default", "grid"),
  methods = "dnn_oracle",
  alpha = 0.05,
  sigma = 0.5,
  clip = 0.001,
  n_folds = 5L,
  verbose = TRUE
)
```

## Arguments

- designs:

  Character vector of DGP names.

- n_vals:

  Integer vector of sample sizes.

- kx_vals:

  Integer vector of covariate dimensions.

- n_reps:

  Number of Monte Carlo replications per cell.

- master_seed:

  Integer seed used to draw the per-replicate seeds. Every cell uses the
  same seed vector so noise realisations are comparable across DGPs,
  sample sizes, and methods (fair comparison).

- k:

  Integer DNN neighbour count. The paper-facing DNN default is `k = 1L`,
  meaning one nearest neighbor is selected inside each size-`s`
  subsample; larger values use the package's generalized top-`k`
  weights. The fixed-`k` extension is essentially direct, but the
  current paper derivations are written for `k = 1`.

- s_fn:

  Function `f(n, kx)` returning the integer kernel order. The
  finite-sample default is \\s = \lfloor 0.4 n^{(k_x+1)/(k_x+2)}
  \rfloor\\.

- tuning:

  Optional kernel-order tuning specification. When supplied, `s` is
  resolved per simulation cell from `tuning_eval(..., mode = "default")`
  and `s_fn` must be omitted.

- tuning_mode:

  One of `"default"` or `"grid"`. `"default"` preserves scalar
  simulation behavior. `"grid"` requires `tuning`; DNN methods emit one
  row per tuning candidate from the integrated
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md)
  grid payload, while non-DNN benchmarks emit one row with
  `candidate_id = NA` and `s = NA`.

- methods:

  Character vector of method names; see
  [`sim_cate_methods()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_methods.md).
  Default: `"dnn_oracle"`.

- alpha:

  Nominal coverage level.

- sigma:

  Noise SD forwarded to the DGP.

- clip:

  Propensity clipping bound forwarded to the DGP and cross-fitting.

- n_folds:

  Number of cross-fitting folds for the `dnn_crossfit` method.

- verbose:

  Print one progress line per cell (default `TRUE`).

## Value

A `data.frame` with raw rows and columns `method`, `design`, `n`, `kx`,
`k`, `tuning_source`, `candidate_id`, `lambda`, `s`, `seed`, `tau_hat`,
`tau0`, `bias`, `sq_err`, `se`, `cover`, `ci_len`, `neff`, localization
and studentization diagnostics described in
[`sim_cate_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_one_rep.md),
and `status`. Failed entries have `NA` for all metric columns.

## Parallelism

When the `future.apply` package is installed, the rep loop within each
cell is parallelised via
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html).
Control the backend by calling
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
before `sim_cate_run()`:

    # All local cores (PC or single SLURM node):
    future::plan(future::multisession)
    raw <- sim_cate_run(...)

Without `future.apply` the function falls back silently to sequential
`lapply`. For multi-node HPC, see
`system.file("paper-replication/cate/hpc", package = "RdistNN")`.

## Fault tolerance

Each replicate is wrapped in `tryCatch`. A total rep failure (e.g. DGP
draw error) records `NA` for all metrics and all methods and emits a
warning; per-method failures are handled inside
[`sim_cate_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_one_rep.md).
[`summarise_sim_results()`](https://JakobJuergens.github.io/RdistNN/reference/summarise_sim_results.md)
reports `n_failed` per candidate-aware group.

## See also

[`sim_cate_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_one_rep.md),
[`summarise_sim_results()`](https://JakobJuergens.github.io/RdistNN/reference/summarise_sim_results.md),
[`sim_cate_methods()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_methods.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single method (default):
raw <- sim_cate_run(designs = "dgp2_nonlinear", n_vals = 500L, n_reps = 10L)

# Multiple methods on the same draws:
future::plan(future::multisession)
raw <- sim_cate_run(
  designs = "dgp2_nonlinear", n_vals = 500L, n_reps = 10L,
  methods = c("dnn_fullsample", "dnn_crossfit", "ols_interact", "causal_forest")
)
future::plan(future::sequential)

summarise_sim_results(raw)
} # }
```
