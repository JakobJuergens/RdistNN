# Print a `summary.plqr_result` (experimental)

Print a `summary.plqr_result` (experimental)

## Usage

``` r
# S3 method for class 'summary.plqr_result'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  A `summary.plqr_result` object.

- digits:

  Number of significant digits for printing.

- ...:

  Not used.

## Value

`x`, invisibly.

## Lifecycle

**Experimental.** PLQR functionality is outside the CATE and regression
workflow focus of the initial public release. Its API and inference
interfaces may change in later releases.

## See also

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).
