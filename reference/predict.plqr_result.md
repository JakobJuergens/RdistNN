# Predict from a `plqr_result` object at new targets

Predict from a `plqr_result` object at new targets

## Usage

``` r
# S3 method for class 'plqr_result'
predict(
  object,
  newdata = NULL,
  X0 = NULL,
  x0 = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `plqr_result` object from
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).

- newdata:

  Optional matrix/data frame of target covariates.

- X0:

  Optional numeric matrix of targets (rows).

- x0:

  Optional numeric vector for a single target.

- se.fit:

  Logical; if `TRUE`, also return direct coefficient standard errors.
  The fitted object must retain density values from `density` or
  `density_col`.

- interval:

  One of `"none"` or `"confidence"`.

- level:

  Confidence level for interval output.

- ...:

  Not used.

## Value

With `interval = "none"` and `se.fit = FALSE`, a coefficient matrix with
one row per PLQR coefficient and one column per target. With
`se.fit = TRUE`, a list with `fit` and `se.fit` matrices. With
`interval = "confidence"`, a list with `fit`, `lwr`, and `upr` matrices,
and also `se.fit` when requested. The coefficient source is attached to
the fit as `attr(., "theta_source")`.

## Details

The returned fit matrix repeats the supplied `theta_eval`
(`meta$theta_source = "supplied"`). Target-specific variation enters
only through the DNN weights and the resulting direct-sandwich standard
errors; no target-specific PLQR root is solved.

Coefficient-level standard errors and intervals require the direct PLQR
sandwich (fit with `density` or `density_col`); a componentwise
moment-level jackknife SE on a fit without a density input is not by
itself a coefficient-level interval — see
`markdown_notes/design/15_plqr_variance_estimation_rationale.md`.
Prediction standard errors always use the analytic sandwich meat,
regardless of the fit's `theta_meat` choice, because the fitted object
does not retain the outcome data needed to re-run the jackknife at new
targets.

## Examples

``` r
set.seed(3)
demo <- plqr_demo_data(n = 80, p = 2, p_d = 2, tau = 0.5, seed = 3)
fit <- plqr_dnn(
  data = demo$data,
  x_cols = c("X1", "X2"),
  D_cols = c("D1", "D2"),
  g_col = "g_oracle",
  mD_cols = c("mD1", "mD2"),
  y_col = "Y",
  theta_eval = demo$truth$theta_true,
  tau = 0.5,
  x0 = c(0, 0),
  method = "exact",
  k = 1L, s = 5L,
  se = FALSE
)
predict(fit, x0 = c(0.2, -0.1))
#>    target_1
#> D1      0.5
#> D2      0.5
#> attr(,"theta_source")
#> [1] "supplied"
```
