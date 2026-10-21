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

Coefficient uncertainty requires a fit with `density` or `density_col`.
`moment_se` describes score variation and is not a coefficient standard
error.

## Lifecycle

**Experimental.** The PLQR API and inference interfaces may change.

## See also

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).
