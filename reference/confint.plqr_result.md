# Confidence intervals for `plqr_result`

Confidence intervals for `plqr_result`

## Usage

``` r
# S3 method for class 'plqr_result'
confint(object, parm = NULL, level = 0.95, target = NULL, ...)
```

## Arguments

- object:

  A `plqr_result` object from
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).

- parm:

  Optional coefficient names or indices.

- level:

  Confidence level.

- target:

  Optional target indices.

- ...:

  Not used.

## Value

A data frame with confidence intervals and a `theta_source` attribute.

## Details

Coefficient-level standard errors and intervals require the direct PLQR
sandwich (fit with `density` or `density_col`); the moment-level
jackknife SE available on fits without a density input has no proven
translation into a coefficient-level interval in this package's source
theory — see
`markdown_notes/design/15_plqr_variance_estimation_rationale.md`.

Intervals are centered at the coefficient value the object represents:
for current fits the supplied `theta_eval`
(`meta$theta_source = "supplied"`), not an estimated local root. The
source is attached as `attr(., "theta_source")`.
