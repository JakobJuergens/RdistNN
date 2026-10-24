# Variance-covariance matrices for `rdistnn_plqr_result` (experimental)

Variance-covariance matrices for `rdistnn_plqr_result` (experimental)

## Usage

``` r
# S3 method for class 'rdistnn_plqr_result'
vcov(object, target = NULL, ...)
```

## Arguments

- object:

  A `rdistnn_plqr_result` object.

- target:

  Optional target index. When omitted, returns one matrix for a
  single-target fit and a named list for a multi-target fit.

- ...:

  Not used.

## Value

A covariance matrix or list of covariance matrices.

## Lifecycle

**Experimental.** The PLQR API and inference interfaces may change.

## See also

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).
