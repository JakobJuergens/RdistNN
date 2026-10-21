# Confidence intervals for `plqr_result` (experimental)

Confidence intervals for `plqr_result` (experimental)

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

Coefficient uncertainty requires a fit with `density` or `density_col`.
`moment_se` describes score variation and is not a coefficient standard
error.

Intervals are centered at the coefficient value the object represents:
for current fits the supplied `theta_eval`
(`meta$theta_source = "supplied"`), not an estimated local root. The
source is attached as `attr(., "theta_source")`.

## Lifecycle

**Experimental.** The PLQR API and inference interfaces may change.

## See also

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).
