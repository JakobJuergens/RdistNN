# Summarise raw IV-GMM simulation results

Aggregates the raw per-replicate output from
[`sim_iv_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_run.md)
to one row per `(method, design, n, kx, q)` cell.

## Usage

``` r
summarise_iv_results(raw)
```

## Arguments

- raw:

  A data frame with columns returned by
  [`sim_iv_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_run.md).

## Value

A data frame with bias, RMSE, mean standard error, coverage, mean CI
length, effective sample size, and overidentification diagnostics.

## See also

[`sim_iv_run()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_run.md)

## Examples

``` r
raw <- sim_iv_run(
  designs = "dgp1_iv_justid_linear",
  n_vals = 100,
  q_vals = c(1, 2),
  n_reps = 2,
  methods = "iv_oracle",
  verbose = FALSE
)
summarise_iv_results(raw)
#>      method                design   n kx q n_reps n_failed       bias      rmse
#> 1 iv_oracle dgp1_iv_justid_linear 100  2 1      2        0 -0.2661551 0.3547985
#>     se_mean cover   ci_len     neff J_median pval_J_median reject_J
#> 1 0.5922282     1 2.321492 15.09131       NA            NA       NA
```
