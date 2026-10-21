# Summarise raw per-replicate CATE simulation results

Aggregates a data.frame of raw replicate-level results (as returned by
[`sim_cate_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_run.md))
to one summary row per
`(method, design, n, kx, candidate_id, lambda, s)` cell and, when
present, per GRF nuisance configuration. Reports bias, RMSE, mean SE,
empirical coverage, mean CI length, empirical sampling dispersion of
`tau_hat`, the ratio of mean SE to empirical dispersion, and
(optionally) mean effective sample size. When the raw rows carry
localization diagnostics, the summary also reports mean realized
localization drift, the mean estimator remainder, physical neighborhood
radius, arm-specific effective sample sizes, and shape diagnostics for
the studentized statistic. Failed or skipped rows are excluded from all
averages and dispersion calculations; their count is reported in
`n_failed`. Here `localized_truth` is the realized DNN-weighted truth,
not the deterministic population branch. Under complete-DNN
exchangeability its expectation equals \\\alpha_s(x_0)\\, so mean
realized localization drift estimates the population branch
displacement. The estimator remainder contains centered score variation
and fitted-state effects.

## Usage

``` r
summarise_sim_results(raw)
```

## Arguments

- raw:

  A `data.frame` with at least the columns `method`, `design`, `n`,
  `kx`, `candidate_id`, `lambda`, `s`, `bias`, `sq_err`, `se`, `cover`,
  `ci_len`, and `status`. Optional columns `neff`, `tau_hat`, and `tau0`
  (returned by
  [`sim_cate_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_run.md))
  are handled gracefully: `neff` is averaged into the summary; `tau_hat`
  is summarized through its empirical standard deviation; `tau0` remains
  a per-rep quantity and is not aggregated. `method_variant`,
  `forest_nuisance`, `y_hat_supplied`, and `w_hat_supplied` extend the
  grouping key. These route columns are required; methods without GRF
  nuisances record typed missing values.

## Value

A `data.frame` with one row per unique group and columns `method`,
`design`, `n`, `kx`, `candidate_id`, `lambda`, `s`, `n_reps`,
`n_failed`, `bias`, `rmse`, `se_mean`, `cover`, `ci_len`, `emp_sd`,
`se_to_sd`, `bias_to_se`, localization, support, and studentized-shape
diagnostics (and `neff` when present in `raw`). The output always
includes the four GRF route fields, with typed missing values when they
do not apply.

## See also

[`sim_cate_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_run.md)

## Examples

``` r
if (FALSE) { # \dontrun{
raw <- sim_cate_run(designs = "dgp1_linear", n_vals = 500L, n_reps = 10L,
                    methods = c("dnn_oracle", "ols_interact"))
summarise_sim_results(raw)
} # }
```
