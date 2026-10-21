# Run the IV-GMM simulation grid

Expands a grid over DGP designs, sample sizes, covariate dimensions, and
instrument counts, then runs
[`sim_iv_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_one_rep.md)
for each replicate. Each cell uses the same replicate seed vector so
methods are compared on aligned draws.

## Usage

``` r
sim_iv_run(
  designs = c("dgp1_iv_justid_linear", "dgp2_iv_homog_overid", "dgp3_iv_heterog_overid",
    "dgp4_iv_weak_local", "dgp5_iv_correlated_instruments",
    "dgp6_iv_nonlinear_firststage"),
  n_vals = c(500L, 1000L, 2000L),
  kx_vals = 2L,
  q_vals = c(1L, 2L, 3L),
  n_reps = 500L,
  master_seed = 42L,
  k = 1L,
  s_fn = function(n, kx) floor(0.4 * n^((kx + 1L)/(kx + 2L))),
  tuning = NULL,
  methods = "dnn_iv_oracle",
  alpha = 0.05,
  sigma = 0.5,
  clip = 0.001,
  n_folds = 5L,
  verbose = TRUE
)
```

## Arguments

- designs:

  Character vector of IV DGP names.

- n_vals:

  Integer vector of sample sizes.

- kx_vals:

  Integer vector of covariate dimensions.

- q_vals:

  Integer vector of requested instrument counts. The just-identified
  design is run once with realized `q = 1`.

- n_reps:

  Number of Monte Carlo replications per cell.

- master_seed:

  Integer seed used locally to generate replicate seeds; the caller's
  saved RNG state is restored. `NULL` uses the current stream.

- k:

  Integer DNN neighbor count.

- s_fn:

  Function `f(n, kx)` returning kernel order. The finite-sample default
  is \\\lfloor 0.4 n^{(k_x+1)/(k_x+2)} \rfloor\\.

- tuning:

  Optional tuning specification. When supplied, omit `s_fn`.

- methods:

  Character vector of methods; see
  [`sim_iv_methods()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_methods.md).

- alpha:

  Nominal test size.

- sigma:

  Noise SD forwarded to
  [`sim_iv_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_data.md).

- clip:

  Clipping bound forwarded to DGPs and nuisance fits.

- n_folds:

  Number of cross-fitting folds.

- verbose:

  Print one progress line per cell.

## Value

Raw replicate-level `data.frame`.

## Parallelism

With `future.apply` installed, replicate execution follows
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
through `future.apply::future_lapply(..., future.seed = TRUE)`.
Otherwise it uses sequential `lapply`.

## Reproducibility

An explicit `master_seed` restores the caller's saved RNG state
(including seed absence) and RNG kind on success or error. `NULL`
consumes the current stream. The hidden Box-Muller normal cache is not
part of R's saved state. Repeating a run requires the same seed, RNG
kind, dependency route, methods, tuning callbacks and execution
configuration.

With `future.apply` installed, `future.seed = TRUE` supplies
L'Ecuyer-CMRG streams to the replicate bodies. Package-controlled draws
agree between sequential and multisession future plans for the same RNG
configuration; external learners and callbacks must also be
reproducible. Replicate-level
[`set.seed()`](https://rdrr.io/r/base/Random.html) inherits the RNG kind
in effect there. Without `future.apply`, replicates instead inherit the
caller's RNG kind, so equal integer seeds need not give equal data draws
across dependency routes. Even with a caller L'Ecuyer-CMRG kind,
stochastic cell-level tuning can differ in later cells because future
mapping and `lapply` advance the caller stream differently. Record
[`RNGkind()`](https://rdrr.io/r/base/Random.html), package versions, the
future plan and optional-dependency availability alongside seeds when
saving a simulation configuration.

## See also

[`sim_iv_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_one_rep.md),
[`summarise_iv_results()`](https://JakobJuergens.github.io/RdistNN/reference/summarise_iv_results.md),
[`sim_iv_methods()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_methods.md)

## Examples

``` r
if (FALSE) { # \dontrun{
raw <- sim_iv_run(
  designs = "dgp2_iv_homog_overid",
  n_vals = 200,
  q_vals = 2,
  n_reps = 10,
  methods = c("dnn_iv_oracle", "iv_oracle")
)
summarise_iv_results(raw)
} # }
```
