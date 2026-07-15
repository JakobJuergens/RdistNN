# Extract the per-target diagnostics frame from a fitted RdistNN object

A tidy wrapper around
[`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md):
the same per-target frame (`n_eff`, weight concentration, distance
radii, tuning provenance, and estimator-specific columns), returned as a
plain `data.frame` with an integer `target_id` consistent with the other
`plot_data_*()` extractors.

## Usage

``` r
plot_data_diagnostics(object)
```

## Arguments

- object:

  A fitted result object with a
  [`diagnostics()`](https://JakobJuergens.github.io/RdistNN/reference/diagnostics.md)
  method.

## Value

A `data.frame` with one row per target point; the original diagnostics
label is kept in `target_label`.

## Examples

``` r
set.seed(1)
n <- 60
x_mat <- matrix(rnorm(n * 2), n, 2)
y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
fit <- np_dnn(
  X = x_mat, x0 = c(0, 0), Y = y_vec,
  k = 1L, s = 10L, se = FALSE
)
plot_data_diagnostics(fit)
#>          target_id  n p k  s method distance    n_eff max_weight
#> target_1         1 60 2 1 10  exact       l2 10.54092  0.1666667
#>          top1_weight_share top5_weight_share top10_weight_share
#> target_1         0.1666667         0.6120561          0.8637521
#>          n_nonzero_weight distance_diagnostics_available radius_rank_k
#> target_1               46                          FALSE            NA
#>          radius_rank_s weighted_mean_distance weighted_median_distance
#> target_1            NA                     NA                       NA
#>          tuning_source tuning_name tuning_default_s tuning_grid_available
#> target_1          <NA>        <NA>               NA                    NA
#>          tuning_grid_n_candidates tuning_lambda_default outcome_available
#> target_1                       NA                    NA              TRUE
#>          weighted_mean_y weighted_sd_y weighted_y_q10 weighted_y_q50
#> target_1     0.002292107     0.4238332     -0.5058686     -0.1309296
#>          weighted_y_q90 target_label
#> target_1       0.586555     target_1
```
