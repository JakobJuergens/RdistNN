# Variance-covariance matrices for `rdistnn_gmm_result`

Variance-covariance matrices for `rdistnn_gmm_result`

## Usage

``` r
# S3 method for class 'rdistnn_gmm_result'
vcov(object, target = NULL, ...)
```

## Arguments

- object:

  A `rdistnn_gmm_result` object.

- target:

  Optional target index. When omitted, returns one matrix for a
  single-target fit and a named list for a multi-target fit.

- ...:

  Not used.

## Value

A covariance matrix or list of covariance matrices.
