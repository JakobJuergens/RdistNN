# Variance-covariance matrices for `plqr_result` (experimental)

Variance-covariance matrices for `plqr_result` (experimental)

## Usage

``` r
# S3 method for class 'plqr_result'
vcov(object, target = NULL, ...)
```

## Arguments

- object:

  A `plqr_result` object.

- target:

  Optional target index. When omitted, returns one matrix for a
  single-target fit and a named list for a multi-target fit.

- ...:

  Not used.

## Value

A covariance matrix or list of covariance matrices.

## Lifecycle

**Experimental.** PLQR functionality is outside the CATE and regression
workflow focus of the initial public release. Its API and inference
interfaces may change in later releases.

## See also

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).
