# Summarize a `plqr_result` (experimental)

Summarize a `plqr_result` (experimental)

## Usage

``` r
# S3 method for class 'plqr_result'
summary(object, level = 0.95, ...)
```

## Arguments

- object:

  A `plqr_result` object from
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).

- level:

  Confidence level.

- ...:

  Not used.

## Value

An object of class `summary.plqr_result`.

## Details

Coefficient-level standard errors and intervals require the direct PLQR
sandwich (fit with `density` or `density_col`); the moment-level
jackknife SE available on fits without a density input has no proven
translation into a coefficient-level interval in this package's source
theory — see
`markdown_notes/design/15_plqr_variance_estimation_rationale.md`.

## Lifecycle

**Experimental.** PLQR functionality is outside the CATE and regression
workflow focus of the initial public release. Its API and inference
interfaces may change in later releases.

## See also

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).
