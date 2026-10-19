# Extract a tidy estimates frame from a fitted RdistNN object

Returns one row per target point and coefficient with the point
estimate, stored jackknife standard error, and symmetric normal
confidence bounds. This is the data layer of the plotting interface:
plot builders consume these frames, and users can consume them directly
for custom figures.

## Usage

``` r
plot_data_estimates(object, level = 0.95, ...)

# Default S3 method
plot_data_estimates(object, level = 0.95, ...)

# S3 method for class 'np_result'
plot_data_estimates(object, level = 0.95, candidates = "default", ...)

# S3 method for class 'cate_result'
plot_data_estimates(object, level = 0.95, candidates = "default", ...)

# S3 method for class 'plqr_result'
plot_data_estimates(object, level = 0.95, quantity = c("moment", "theta"), ...)

# S3 method for class 'gmm_result'
plot_data_estimates(object, level = 0.95, ...)

# S3 method for class 'affine_ratio_result'
plot_data_estimates(object, level = 0.95, candidates = "default", ...)
```

## Arguments

- object:

  A fitted result object (`np_result`, `cate_result`,
  `affine_ratio_result`, `clate_result`, `plqr_result`, or `gmm_result`;
  `tdnn_np_result` dispatches through `np_result`).

- level:

  Confidence level for the normal bounds (default `0.95`).

- ...:

  Passed to methods; not used.

- candidates:

  For grid-aware methods, `"default"`, `"all"`, or a character vector of
  stored candidate ids.

- quantity:

  For `plqr_result` fits, which quantity to extract: `"moment"`
  (default) uses the localized PLQR moments and their jackknife
  `moment_se`; `"theta"` uses the supplied coefficient values (repeated
  `theta_eval`) and the direct-sandwich `theta_se`. While no root solver
  exists, moments are the target-specific quantity, so they are the
  default. The choice is recorded in `attr(., "quantity")`.

## Value

A `data.frame` with columns `target_id` (integer, `1..m`), `coef`
(character; `"mu"`/`"tau"`/`"theta"` for scalar fits, moment or
coefficient names for vector fits), `estimate`, `se`, `ci_lower`,
`ci_upper`, and — when the fit stores evaluation points — target
coordinates `x0_1..x0_p`, with the original covariate names recorded in
`attr(., "x0_names")`. Fits without stored standard errors return `NA`
in `se` and the bound columns. The confidence level is recorded in
`attr(., "level")`. Ratio fits also retain unavailable Wald bounds as
`NA` when the ratio is undefined or its score variance is zero.

## Experimental methods

Methods for `plqr_result` objects are experimental, as is all PLQR
functionality in this release; see
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).

## Examples

``` r
set.seed(1)
n <- 60
x_mat <- matrix(rnorm(n * 2), n, 2)
y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
fit <- np_dnn(
  X = x_mat, X0 = rbind(c(0, 0), c(1, 0)), Y = y_vec,
  k = 1L, s = 10L, se = FALSE
)
plot_data_estimates(fit)
#>   target_id coef    estimate se ci_lower ci_upper x0_1 x0_2
#> 1         1   mu 0.002292107 NA       NA       NA    0    0
#> 2         2   mu 0.832118414 NA       NA       NA    1    0
```
