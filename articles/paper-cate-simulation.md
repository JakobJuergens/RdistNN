# CATE Simulation Study

## Overview

This vignette documents the simulation infrastructure for the companion
conditional Z-estimation paper listed by `citation("RdistNN")`. It shows
how to use the package-level simulation functions on a small synthetic
example. The full Monte Carlo run for the paper main grid is in
`inst/paper-replication/cate/run.R`.

The manuscript repo `Cond_Z_Estimation/Codes/Simulations/` now provides
the paper-facing orchestration layer on top of these package functions:

- chunked resumable runs with run manifests and generated run READMEs,
- the separate nuisance-efficiency launcher, and
- fixed-`\lambda` sweep wrappers that write one output folder per
  `\lambda`.

To access that script after installing the package:

``` r

system.file("paper-replication/cate/run.R", package = "RdistNN")
```

## Data-generating processes

Six DGPs are available via
[`sim_cate_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_data.md):

| Design | `tau(x)` | Propensity | Covariate design | Notes |
|----|----|----|----|----|
| `dgp1_linear` | `0.5 x1 - 0.3 x2` | logistic | iid N(0,I) | correctly specified for OLS-interaction |
| `dgp2_nonlinear` | `0.5 + 0.75 sin(1.5 x1)` | logistic | iid N(0,I) | sharper local misspecification design |
| `dgp3_sparse_overlap` | same as DGP 2 | near-boundary except in a tight pocket at `x0` | iid N(0,I) | stress-tests local ESS adaptation |
| `dgp4_clustered` | same as DGP 2 | logistic | mixture: 80% N(x0, 0.2²I) + 10% N((4.5,2), 1.25²I) + 10% N((-2.5,-2), 1.25²I) | density-adaptivity DGP; DNN shrinks bandwidth at dense `x0` |
| `dgp5_correlated` | `sin(2x1) + 0.5 x1²` | steeper logistic | N(0,Σ), ρ=0.85; x0=(0.75,0) | shape-adaptivity DGP; `dnn_oracle_mahal` corrects for off-diagonal covariance |
| `dgp6_nuisance_efficiency` | `0.5 + 0.60 sin(1.5 x1)` | logistic, simple benchmark nuisance structure | iid N(0,I) | focused `dnn_fullsample` versus `dnn_crossfit` design |

All DGPs return oracle nuisances (`e_oracle`, `m0_oracle`, `m1_oracle`,
`tau_oracle`), the paper-facing oracle AIPW score (`aipw_oracle`), and
the legacy Robinson pseudo-outcome (`rho_oracle`).

``` r

dgp <- sim_cate_data("dgp2_nonlinear", n = 300, kx = 2, seed = 42)
head(dgp$data[, c("X1", "X2", "D", "Y", "tau_oracle", "e_oracle")])
#>           X1           X2 D          Y  tau_oracle  e_oracle
#> 1  1.3709584 -0.004620768 1  2.7348310  1.16328183 0.6654750
#> 2 -0.5646982  0.760242168 0 -1.2998322 -0.06199628 0.3401849
#> 3  0.3631284  0.038990913 0  1.1537875  0.88861641 0.5404287
#> 4  0.6328626  0.735072142 0  0.9851532  1.10975343 0.4872266
#> 5  0.4042683 -0.146472627 0  0.4118522  0.92743630 0.5684109
#> 6 -0.1061245 -0.057887335 1  1.4781386  0.38111351 0.4939706
dgp$tau_x0   # true CATE at the default target x0 = (1, 0)
#> [1] 1.248121
```

## Single replicate

[`sim_cate_one_rep()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_one_rep.md)
draws one dataset, fits the DNN Z-estimator using the AIPW score, and
returns per-replicate metrics:

``` r

tuning_default <- tuning_dnn_paper()
s_default <- tuning_eval(tuning_default, n = 300, kx = 2, k = 1)$s

r <- sim_cate_one_rep(
  n      = 300,
  kx     = 2,
  k      = 1,
  s      = s_default,
  design = "dgp2_nonlinear",
  seed   = 1
)
r
#>       method candidate_id lambda  s  tau_hat     tau0      bias      sq_err
#> 1 dnn_oracle         <NA>     NA 28 1.289072 1.248121 0.0409505 0.001676943
#>          se cover    ci_len     neff localized_truth branch_drift
#> 1 0.1805794     1 0.7078582 20.09801        1.196847  -0.05127403
#>   estimator_remainder studentized_error score_sd0 standardized_branch_drift
#> 1          0.09222453         0.2267728  1.031413                -0.2228647
#>   physical_radius_90 treated_ess control_ess top5_score_share status
#> 1          0.5221937    12.00479     8.38294         0.457057     ok
```

## Small Monte Carlo run

The following runs 5 replications per cell over a tiny grid to
illustrate the workflow. The manuscript-facing main production grid uses
`n_reps = 10000`, `n_vals = c(100, 250, 500, 1000)`, and the retained
designs `dgp1_linear`, `dgp2_nonlinear`, `dgp4_clustered`, and
`dgp5_correlated`. The focused nuisance-efficiency comparison and the
fixed-`\lambda` sweeps are driven from the manuscript repo rather than
this package vignette.

``` r

raw <- sim_cate_run(
  designs     = c("dgp1_linear", "dgp2_nonlinear"),
  n_vals      = 300L,
  kx_vals     = 2L,
  n_reps      = 5L,
  master_seed = 42L,
  tuning      = tuning_dnn_paper(),
  verbose     = FALSE
)
head(raw)
#>        method         design   n kx k tuning_source candidate_id lambda  s
#> 1  dnn_oracle    dgp1_linear 300  2 1        tuning         <NA>    0.4 28
#> 11 dnn_oracle    dgp1_linear 300  2 1        tuning         <NA>    0.4 28
#> 12 dnn_oracle    dgp1_linear 300  2 1        tuning         <NA>    0.4 28
#> 13 dnn_oracle    dgp1_linear 300  2 1        tuning         <NA>    0.4 28
#> 14 dnn_oracle    dgp1_linear 300  2 1        tuning         <NA>    0.4 28
#> 15 dnn_oracle dgp2_nonlinear 300  2 1        tuning         <NA>    0.4 28
#>      seed     tau_hat     tau0         bias       sq_err        se cover
#> 1   61413  0.30162522 0.000000  0.301625221 9.097777e-02 0.1814035     1
#> 11  54425  0.21851858 0.000000  0.218518582 4.775037e-02 0.1614648     1
#> 12 623844 -0.02943035 0.000000 -0.029430347 8.661454e-04 0.1537416     1
#> 13  74362 -0.21014204 0.000000 -0.210142045 4.415968e-02 0.2143997     1
#> 14  46208 -0.26902641 0.000000 -0.269026414 7.237521e-02 0.2089030     1
#> 15  61413  1.23930527 1.248121 -0.008815971 7.772135e-05 0.1525520     1
#>       ci_len     neff localized_truth branch_drift estimator_remainder
#> 1  0.7110885 20.09801     0.003754586  0.003754586          0.29787064
#> 11 0.6329303 20.09801     0.012698141  0.012698141          0.20582044
#> 12 0.6026559 20.09801    -0.005642472 -0.005642472         -0.02378788
#> 13 0.8404312 20.09801    -0.002948573 -0.002948573         -0.20719347
#> 14 0.8188848 20.09801     0.023634586  0.023634586         -0.29266100
#> 15 0.5979929 20.09801     1.156890838 -0.091230402          0.08241443
#>    studentized_error score_sd0 standardized_branch_drift physical_radius_90
#> 1         1.66273129  1.000000                0.01683211          0.3477558
#> 11        1.35335137  1.000000                0.05692679          0.4306427
#> 12       -0.19142739  1.000000               -0.02529566          0.4041881
#> 13       -0.98014170  1.000000               -0.01321869          0.4067113
#> 14       -1.28780521  1.000000                0.10595575          0.3700703
#> 15       -0.05778994  1.031413               -0.39653678          0.5580878
#>    treated_ess control_ess top5_score_share status
#> 1     8.425537   11.691286        0.4454872     ok
#> 11   10.823350    9.396065        0.4710436     ok
#> 12    7.820847   12.277236        0.3730882     ok
#> 13   13.357003    7.178954        0.5301387     ok
#> 14   12.869978    9.113387        0.4887371     ok
#> 15   13.286536    7.069731        0.4273958     ok
```

[`summarise_sim_results()`](https://JakobJuergens.github.io/RdistNN/reference/summarise_sim_results.md)
aggregates raw per-replicate rows to one summary row per cell:

``` r

tab <- summarise_sim_results(raw)
print(tab, digits = 3)
#>       method         design   n kx candidate_id lambda  s n_reps n_failed
#> 1 dnn_oracle    dgp1_linear 300  2         <NA>    0.4 28      5        0
#> 2 dnn_oracle dgp2_nonlinear 300  2         <NA>    0.4 28      5        0
#>       bias  rmse se_mean cover ci_len emp_sd se_to_sd bias_to_se
#> 1  0.00231 0.226   0.184     1  0.721  0.253    0.727     0.0126
#> 2 -0.03949 0.206   0.209     1  0.821  0.226    0.926    -0.1885
#>   localized_truth_mean branch_drift estimator_remainder score_sd0
#> 1               0.0063       0.0063            -0.00399      1.00
#> 2               1.1899      -0.0582             0.01871      1.03
#>   standardized_branch_drift physical_radius_90_mean physical_radius_90_p25
#> 1                    0.0282                   0.392                  0.363
#> 2                   -0.2530                   0.496                  0.437
#>   physical_radius_90_p50 physical_radius_90_p75 treated_ess control_ess
#> 1                  0.404                  0.415        10.7        9.93
#> 2                  0.465                  0.567        11.6        9.00
#>   top5_score_share studentized_mean studentized_sd studentized_skewness
#> 1            0.462            0.111          1.341                0.184
#> 2            0.522           -0.200          0.928               -0.759
#>   studentized_excess_kurtosis studentized_q025 studentized_q05 studentized_q50
#> 1                      -1.701            -1.29           -1.29         -0.1914
#> 2                      -0.358            -1.70           -1.70         -0.0578
#>   studentized_q95 studentized_q975 studentized_abs_q95 neff
#> 1           1.663            1.663                1.66 20.1
#> 2           0.853            0.853                1.70 20.1
```

## Parallel execution

[`sim_cate_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_run.md)
uses
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
for the rep loop when `future.apply` is installed, and falls back
silently to sequential `lapply` otherwise. The execution backend is
controlled by calling
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
before the run — no changes to
[`sim_cate_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_run.md)
are needed.

### Personal computer

``` r

library(future)
library(future.apply)

plan(multisession)            # use all local cores
raw <- sim_cate_run(
  designs = "dgp2_nonlinear",
  n_vals  = 500L,
  n_reps  = 50L
)
plan(sequential)              # reset when done
```

### SLURM HPC (array jobs)

For the full paper main grid on a cluster, use the HPC scripts in
`inst/paper-replication/cate/hpc/`. Each SLURM array task runs one
`(design, n, kx)` cell (all 10,000 reps), and uses the cores allocated
by `--cpus-per-task` for within-task rep parallelism.

``` r

# Locate the HPC scripts after installing the package:
system.file("paper-replication/cate/hpc", package = "RdistNN")
```

Workflow on the cluster:

``` bash
# 1. Install RdistNN on the cluster, then from your working directory:
sbatch /path/to/hpc/submit.sh        # submits 16 array tasks (edit --array if grid changes)

# 2. After all tasks complete, merge the per-cell results:
Rscript /path/to/hpc/merge.R
```

Results are written to `./results/` by default. Override with
`export SIM_OUTDIR=/your/path` before calling `sbatch` or `merge.R`.

## Full paper replication (single machine)

To run the full simulation sequentially (or with `plan(multisession)`)
on one machine:

``` r

source(system.file("paper-replication/cate/run.R", package = "RdistNN"))
```

Results are written to `./results/sim_cate_raw.rds` and
`./results/sim_cate_table.csv` in the current working directory.

This package-side script is the manuscript-facing **main-grid** baseline
only. If you want the chunked manuscript output with manifests, the
nuisance-efficiency run, or the `\lambda`-sweep folders, use the
launchers in `Cond_Z_Estimation/Codes/Simulations/`.

## Competitor benchmarks

Stub implementations of causal forest, local linear forest, naive OLS,
and interaction OLS are in `inst/paper-replication/cate/benchmarks.R`.
These are not part of the package core; they use
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) and require
optional packages (`grf`, etc.) listed under `Suggests` in
`DESCRIPTION`.

``` r

system.file("paper-replication/cate/benchmarks.R", package = "RdistNN")
```
