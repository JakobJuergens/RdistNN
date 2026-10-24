# Extract coefficients from a `rdistnn_plqr_result` (experimental)

Returns the coefficient value represented by the object. For current
[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md)
fits this is the supplied `theta_eval` repeated across targets
(`meta$theta_source = "supplied"`); no local PLQR root is solved. The
source is attached as `attr(., "theta_source")`.

## Usage

``` r
# S3 method for class 'rdistnn_plqr_result'
coef(object, ...)
```

## Arguments

- object:

  A `rdistnn_plqr_result` object.

- ...:

  Not used.

## Value

Matrix of PLQR coefficient values with a `theta_source` attribute.

## Lifecycle

**Experimental.** The PLQR API and inference interfaces may change.

## See also

[`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).
