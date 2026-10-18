# Compute diagnostics for fitted localized estimators

Returns a compact diagnostics table for supported fit objects.
Per-target DNN diagnostics report fixed top-1, top-5, and top-10 weight
concentration shares. Distance diagnostics such as `radius_rank_k` and
`radius_rank_s` are ordinary nearest-rank scale summaries;
`radius_rank_s` is not the DNN support radius. For exact DNN fits the
finite-population support edge is rank `n - s + k`.

## Usage

``` r
# S3 method for class 'cate_result'
diagnostics(object, candidates = "default", ...)

diagnostics(object, candidates = "default", ...)

# Default S3 method
diagnostics(object, candidates = "default", ...)

# S3 method for class 'dnn_diagnostics'
print(x, ...)

# S3 method for class 'np_result'
diagnostics(object, candidates = "default", ...)

# S3 method for class 'plqr_result'
diagnostics(object, candidates = "default", ...)
```

## Arguments

- object:

  A fitted estimator object, currently an object returned by
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md),
  [`affine_ratio_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/affine_ratio_dnn.md)
  or
  [`clate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/clate_dnn.md).

- candidates:

  Candidate rows to report for grid fits. The default `"default"`
  preserves the per-target diagnostics table for the selected fit. Use
  `"all"` or a character vector of candidate ids to report
  candidate-level grid diagnostics from `object$tuning_grid`. Candidate
  grid diagnostics include `top5_weight_share`; the cutoff is recorded
  as `top_n = 5` and is independent of the estimator's `k`.

- ...:

  Not used.

- x:

  A `dnn_diagnostics` object.

## Value

An object of class `dnn_diagnostics`.

`x`, invisibly.

## Details

PLQR fits additionally report residualized-design Gram conditioning
(`gram_min_eig`, `gram_condition_number`) and, when direct-sandwich
coefficient inference is available, target-level Jacobian conditioning:
`jacobian_available`, `jacobian_min_singular_value`,
`jacobian_condition_number`, `jacobian_rank`, density summaries over the
weight support (`density_min_support`, `density_max_support`, plus the
weighted `density_mean_weighted` and `density_median_weighted`), and
`theta_vcov_available`. A target can have an acceptable residualized
Gram but an ill-conditioned density-weighted Jacobian, so
coefficient-inference reliability should be judged from the Jacobian
block.

## Experimental methods

Methods for `plqr_result` objects are experimental, as is all PLQR
functionality in this release; see
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).

## Examples

``` r
set.seed(1)
x_mat <- matrix(rnorm(40 * 2), 40, 2)
y_vec <- x_mat[, 1] + rnorm(40, sd = 0.2)
fit_np <- np_dnn(X = x_mat, x0 = c(0, 0), Y = y_vec, se = FALSE)
diagnostics(fit_np)
#> <dnn_diagnostics>
#>   targets: 1
#>   n_eff  : min=11.4527, median=11.4527, max=11.4527
#>   weights: max(top5 share)=0.5771, max(weight)=0.15
#>   dist   : available for 0/1 target(s)
#>   outcome: available for 1/1 target(s)

demo_cate <- cate_demo_data(n = 50, p = 2, seed = 2)
fit_cate <- cate_dnn(
  data = demo_cate$data,
  x_cols = c("X1", "X2"),
  y_col = "Y",
  t_col = "D",
  m_col = "m_oracle",
  e_col = "e_oracle",
  x0 = demo_cate$x0,
  se = FALSE
)
diagnostics(fit_cate)
#> <dnn_diagnostics>
#>   targets: 1
#>   n_eff  : min=12.4648, median=12.4648, max=12.4648
#>   weights: max(top5 share)=0.5457, max(weight)=0.14
#>   dist   : available for 0/1 target(s)
#>   overlap: available for 1/1 target(s)

demo_plqr <- plqr_demo_data(n = 60, p = 2, p_d = 1, seed = 3)
fit_plqr <- plqr_dnn(
  data = demo_plqr$data,
  x_cols = c("X1", "X2"),
  D_cols = "D1",
  g_col = "g_oracle",
  mD_cols = "mD1",
  y_col = "Y",
  theta_eval = demo_plqr$truth$theta_true,
  x0 = demo_plqr$x0,
  se = FALSE
)
diagnostics(fit_plqr)
#> <dnn_diagnostics>
#>   targets: 1
#>   n_eff  : min=13.2383, median=13.2383, max=13.2383
#>   weights: max(top5 share)=0.5241, max(weight)=0.1333
#>   dist   : available for 0/1 target(s)
#>   gram   : available for 1/1 target(s)
```
