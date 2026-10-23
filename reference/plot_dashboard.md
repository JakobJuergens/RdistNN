# Compose a diagnostics dashboard for a fitted object

A patchwork composition of the estimate panel with the weight and
overlap diagnostics, annotated with the fit's metadata (`n`, `m`, `k`,
`s` with provenance, distance, and jackknife type). CATE and CLATE
dashboards use estimates, overlap mass, local propensity, and cumulative
weight share; CLATE overlap describes the instrument. Stored tuning
grids add a sensitivity panel. nonparametric-regression dashboards use
estimates, the diagnostics facets, and the weight decay and
concentration views. Panels whose optional retained inputs are absent
are substituted (the diagnostics facets take the place of a missing
overlap panel) or dropped. Invalid plot options, malformed retained
state and unexpected rendering failures propagate.

## Usage

``` r
plot_dashboard(object, ...)

# Default S3 method
plot_dashboard(object, ...)

# S3 method for class 'cate_result'
plot_dashboard(object, ...)

# S3 method for class 'clate_result'
plot_dashboard(object, ...)

# S3 method for class 'np_result'
plot_dashboard(object, ...)
```

## Arguments

- object:

  A fitted result object (`cate_result`, `clate_result` or `np_result`).

- ...:

  Passed on to
  [`plot_estimates()`](https://JakobJuergens.github.io/RdistNN/reference/plot_estimates.md)
  for the estimates panel.

## Value

A patchwork object (or a named list of ggplots when patchwork is
unavailable).

## Details

When patchwork is not installed the panels are returned as a named list
of ggplot objects with a message.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)) {
  set.seed(1)
  n <- 100
  x_mat <- matrix(rnorm(n * 2), n, 2)
  y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
  fit <- np_dnn(
    X = x_mat, X0 = cbind(seq(-1, 1, length.out = 6), 0), Y = y_vec,
    k = 1L, s = 12L, return = c("W", "dx")
  )
  plot_dashboard(fit)
}
```
