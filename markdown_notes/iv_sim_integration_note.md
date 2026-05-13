# IV simulation integration for `RdistNN`

Date: 2026-05-13

## Bottom line

The clean way to add the IV-type example simulations from
`Cond_GMM_with_DNN` to `RdistNN` is to mirror the existing
`sim_cate_*` simulation layer one-to-one, *not* to bolt IV-specific
logic onto `sim_cate_*` or to inline it next to the estimator.

This note covers only the **simulation and demo-data layer**. The
estimator architecture — `gmm_dnn()`, the moment-model object, the
identity/efficient two-step, sandwich SEs and the local `J`-test —
is the subject of `local_gmm_rdistnn_implementation_note.md`
(2026-03-30) and is not retread here. Anything that depends on the
estimator landing is flagged below as **gmm_dnn-gated**.

The target paper section is
`Cond_GMM_with_DNN/Paper/Parts/05_Applications/` (IV-CATE template,
IV simulation, empirical IV slot). The DR-IV moment, local Jacobian,
and local `J`-statistic are pinned down in
`Cond_GMM_with_DNN/Paper/Parts/02_Setup_and_Conditions/02a_Setup_CondGMM_Problem.tex`,
example `ex:iv_cate`.

## Scope split

The work decomposes cleanly into two phases. Phase A can land now in
the same PR or before `gmm_dnn()`; phase B is the integration step
once the estimator exists.

### Phase A — can land now

These files are estimator-free. They take ordinary
`(Y, D, Z, X)` data and produce ground-truth quantities or
oracle/2SLS comparisons.

- `R/sim-iv-data.R` — DGP dispatcher and the six IV DGPs (catalog
  below). Returns a list with `data`, `x0`, the local target
  (`tau_x0` or `theta_x0`, naming flagged below), `q`, and a
  `truth` slot with oracle nuisances.
- `R/demo-data-iv.R` — small `iv_demo_data()` and `iv_demo_sets()`
  factories paralleling `cate_demo_data()` / `cate_demo_sets()`,
  with `inject` / `inject_col` / `idx_inject` hooks for non-finite
  injection testing of downstream IV code.
- Oracle / classical comparators in `R/sim-iv-benchmarks.R`:
  - `.method_iv_oracle` — uses oracle `(π_1,…,π_q, μ_Y, μ_D)` and
    weights, computes the analytic local efficient `τ̂(x_0)` and
    oracle local `J`.
  - `.method_2sls_global` — global 2SLS; serves as a sanity floor.
  - `.method_2sls_interact` — interaction 2SLS in `x` (the IV
    analogue of `ols_interact` in the CATE grid).
  - `.method_local_2sls` — fixed-bandwidth local 2SLS, mirroring
    `.method_local_constant`.
- Tests in `tests/testthat/test-sim-iv-data.R` and
  `test-sim-iv-benchmarks.R` mirroring `test-sim-cate*.R`.

### Phase B — gmm_dnn-gated

- `R/sim-iv-run.R` — `sim_iv_one_rep()` and `sim_iv_run()`. Grid is
  `design × n × kx × q × n_reps`. `master_seed` produces an aligned
  seed vector so methods see the same draws. Same `tryCatch`-with-
  NA-rows failure handling as `sim_cate_run()`.
- Method-registry entries `.method_dnn_iv_fullsample`,
  `.method_dnn_iv_crossfit`, and (optionally)
  `.method_dnn_iv_oracle_mahal`, wiring through `gmm_dnn()` with
  the DR-IV moment object.
- `R/sim-iv-summary.R` — `summarise_iv_results(raw)` aggregating to
  `(method, design, n, kx, q)` cells. Adds the four `J`-test
  columns (see "Metric extension" below) on top of the
  `sim_cate_*` schema.
- `inst/paper-replication/iv/run.R` and `hpc/` mirroring
  `inst/paper-replication/cate/`.

## File layout

```
R/sim-iv-data.R          # DGP catalog, oracle nuisances
R/sim-iv-run.R           # sim_iv_one_rep(), sim_iv_run()           [gated]
R/sim-iv-benchmarks.R    # .method_* family for IV
R/sim-iv-summary.R       # summarise_iv_results()                   [gated]
R/demo-data-iv.R         # iv_demo_data(), iv_demo_sets()
inst/paper-replication/iv/run.R    # PAPER_GRID, paper-facing run   [gated]
inst/paper-replication/iv/hpc/     # submission wrappers            [gated]
tests/testthat/test-sim-iv-data.R
tests/testthat/test-sim-iv-benchmarks.R
tests/testthat/test-sim-iv.R                                        [gated]
tests/testthat/test-demo-data-iv.R
```

Naming chosen so `grep ^sim_iv` and `grep ^iv_demo_` partition the
new surface cleanly and stay disjoint from `sim_cate_*` /
`cate_demo_*`.

## `.sim_metric()` extension

`sim_cate_*` currently builds each per-rep row through
`.sim_metric(tau_hat, se_hat, tau0, alpha, neff)`, which fixes the
columns `tau_hat, tau0, bias, sq_err, se, cover, ci_len, neff`.

For IV with possibly `q > p` we need overidentification diagnostics
first-class. The cleanest path is a new helper
`.sim_metric_iv(theta_hat, se_hat, theta0, alpha, neff,
                J = NA_real_, df = NA_integer_)`
that returns the existing CATE bundle plus

- `J`        — realized local J-statistic at `x_0`
- `df`       — `q - p` (integer; `NA` when overidentification is
                vacuous, e.g. `q = p = 1`)
- `pval_J`   — `pchisq(J, df, lower.tail = FALSE)` with `df > 0`
                guard
- `reject_J` — `as.integer(pval_J < alpha)` with `NA` propagation

Reusing the same scalar-target columns (`theta_hat, theta0, bias,
se, cover, ci_len`) keeps `summarise_iv_results()` mechanically
similar to `summarise_sim_results()`; the J columns aggregate by
`mean(reject_J)` (size/power) and `median(J)` / `median(pval_J)`
(diagnostics).

For the oracle and 2SLS benchmarks that lack a `J` we just pass
`J = NA_real_, df = NA_integer_`; this keeps a uniform row schema
across all methods.

## DGP catalog (six designs)

The catalog mirrors the size of the CATE catalog but reorganizes
around overidentification, which is the structural feature this
paper studies. Each design fixes one phenomenon and varies it cleanly
with `n`, `kx`, and `q`.

1. **`dgp1_iv_justid_linear`** — `q = 1`, linear first stage and
   reduced form. Wald-IV sanity floor; oracle, global 2SLS,
   interaction 2SLS, local 2SLS, and DNN should all coincide
   asymptotically. `df = 0`, J-columns are `NA`.
2. **`dgp2_iv_homog_overid`** — `q ∈ {2, 3}` valid instruments, all
   identifying the same `τ(x_0)`. Size study for the local `J`-test.
3. **`dgp3_iv_heterog_overid`** — `q ∈ {2, 3}` with deliberately
   different complier LATEs at `x_0`. Power study for the local
   `J`-test; should also expose efficiency loss of the equally
   weighted (identity-first-step) estimator versus the efficient
   second step.
4. **`dgp4_iv_weak_local`** — first-stage strength `π_k(x)` shrinks
   to zero near `x_0`. Diagnostic against local many-weak-IV failure
   modes; expected to stress sandwich SE coverage and J-test size.
5. **`dgp5_iv_correlated_instruments`** — overidentified with
   highly correlated instruments (`Corr(Z_j, Z_k)` near 1). Tests
   stability of `Ω̂(x)^{-1}` and the efficient-step inversion.
6. **`dgp6_iv_nonlinear_firststage`** — nonlinear `π_k(x)` and
   `μ_Y(x)`, `μ_D(x)`. The headline DNN-versus-classical comparison;
   parallels `dgp2_nonlinear` in the CATE grid.

Each design returns oracle nuisances under the column convention

```
piZ1_oracle, ..., piZq_oracle,   # Pr(Z_k = 1 | X) or E[Z_k | X]
muY_oracle,                       # E[Y | X]
muD_oracle,                       # E[D | X]
tau_oracle (or theta_oracle)      # τ(X) at every row
```

so the oracle benchmark can construct
`g(W; τ, η) = (Z_k - π_k)(Y - μ_Y - (D - μ_D)τ)` directly without
re-deriving the DR form per DGP. Binary instruments are encouraged
in designs 1–4; design 5 uses continuous correlated `Z`; design 6
uses one binary and the rest continuous to exercise mixed types.

## `q` as a new grid axis

`sim_cate_run()` iterates over `design × n × kx`. `sim_iv_run()`
adds `q` as an outermost grid axis. The natural sweep is
`q_vals = c(1L, 2L, 3L)` with most rows at `q = 2`. Designs that are
just-identified by construction (`dgp1_iv_justid_linear`) ignore
their `q` argument and emit a single row at `q = 1`; the run loop
must skip duplicate `q` cells for such designs rather than running
them three times.

## Paper-replication grid

`inst/paper-replication/iv/run.R` should hold a `PAPER_GRID` list
matching the style of `inst/paper-replication/cate/run.R`. A
reasonable phase-B starting grid (refine when the paper figure
tables are drafted):

```r
PAPER_GRID <- list(
  designs     = c("dgp1_iv_justid_linear",
                  "dgp2_iv_homog_overid",
                  "dgp3_iv_heterog_overid",
                  "dgp5_iv_correlated_instruments",
                  "dgp6_iv_nonlinear_firststage"),
  n_vals      = c(250L, 500L, 1000L, 2000L),
  kx_vals     = 2L,
  q_vals      = c(2L, 3L),
  n_reps      = 10000L,
  master_seed = 42L,
  k           = 1L,
  s_fn        = function(n, kx) floor(n ^ ((kx + 1L) / (kx + 2L))),
  methods     = c("dnn_iv_fullsample", "dnn_iv_crossfit",
                  "local_2sls", "2sls_interact"),
  alpha       = 0.05,
  clip        = 1e-3,
  n_folds     = 5L
)
```

`dgp4_iv_weak_local` is intentionally a diagnostic side-table, not
part of the main grid, paralleling how
`dgp6_nuisance_efficiency` is treated in the CATE replication.

## Testing strategy

Mirroring `test-sim-cate*.R`:

- `test-sim-iv-data.R` — for each DGP: returned object shape,
  presence of oracle nuisance columns, no `NA` in oracle slots,
  `dim(data)` matches `n`, `tau_x0` (or `theta_x0`) finite,
  binary-instrument designs have `Z ∈ {0, 1}`.
- `test-sim-iv-benchmarks.R` — for each `.method_iv_*`:
  deterministic given seed, returns the full
  `.sim_metric_iv()` schema, J columns sensible
  (`NA` for `q = 1`, finite and `df = q - 1` otherwise).
- `test-demo-data-iv.R` — `iv_demo_data()` is dual-interface
  (matrix / data.frame), injection hooks honored, `iv_demo_sets()`
  returns the expected named fixtures.
- (gated) `test-sim-iv.R` — `sim_iv_one_rep()` and a tiny
  `sim_iv_run()` smoke test at `n_reps = 4`, `master_seed = 1L`,
  asserts the row schema and that rows are aligned across methods
  by seed.

Determinism follows the existing pattern: per-rep seeds are derived
from `master_seed` via the same helper used by `sim_cate_run()`,
not from `set.seed()` at the top of the script.

## Open question — local-IV target naming

The Cond_GMM paper uses `τ(x)` in `ex:iv_cate` (a local LATE), but
the wider local-GMM machinery uses `α(x)` or `θ(x)` as the generic
local parameter. The simulation layer should commit to **one** name
once `gmm_dnn()` lands, and pick the same name as the estimator
return slot. Until then, `sim-iv-data.R` should expose both
`theta_x0` and a `tau_x0 = theta_x0` alias and let downstream code
decide; that lets us write the DGPs now without a forward-rename
later. Flag, do not commit, the choice in this note.

## Vignette

A short `vignettes/local-iv-estimation.Rmd` parallels
`vignettes/conditional-z-estimation.Rmd`. Gated on `gmm_dnn()` so
the vignette can actually compute the local DR-IV estimator. Until
then, keep the IV exposition inside the estimator design note, not
in vignettes — design note 11 (examples/vignettes/pkgdown policy)
specifically wants vignettes to be runnable on a clean install.

## Summary of integration steps

1. Land Phase A: DGPs, demo data, oracle/2SLS benchmarks, tests.
   This is independent of `gmm_dnn()` and gives us a working
   simulation harness for the classical comparators immediately.
2. Land `gmm_dnn()` per `local_gmm_rdistnn_implementation_note.md`.
3. Land Phase B: register `.method_dnn_iv_*`, wire
   `sim_iv_one_rep()`, write `summarise_iv_results()`, and add
   `inst/paper-replication/iv/`.
4. Resolve the `τ` vs `θ` naming question against the final
   `gmm_dnn()` return slot.
5. Cut the IV vignette.

## Implementation status as of 2026-05-13

The package now contains the full first IV simulation stack described above:

- Phase A is implemented: IV DGPs, demo-data helpers, oracle/classical
  benchmarks, and tests.
- The `gmm_dnn()`-gated estimator layer is implemented: `new_gmm_model()`,
  `gmm_model_dr_iv()`, `gmm_dnn()`, and S3 methods for `gmm_result`.
- Phase B is implemented: `sim_iv_one_rep()`, `sim_iv_run()`, DNN IV method
  registry entries, and `summarise_iv_results()`.
- The paper-facing replication layer is implemented under
  `inst/paper-replication/iv/`, with a shared grid helper, single-machine run
  script, SLURM task/merge scripts, and a smoke-mode grid for package tests.
- The runnable local IV article is implemented as
  `vignettes/local-iv-estimation.Rmd`.

Remaining IV simulation work is now paper-output refinement rather than package
surface construction: add any diagnostic weak-IV side tables, tune cluster
resources after timing production runs, and update manuscript-facing table or
figure code once the preferred output layout is fixed.
