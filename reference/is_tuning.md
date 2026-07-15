# Is this a tuning spec?

A tuning spec is a list with an `eval(n, kx, k, se, jackknife, ctx)`
function and optional `eval_grid(...)` for candidate generation.

## Usage

``` r
is_tuning(x)
```

## Arguments

- x:

  Any R object.

## Value

Logical.

## Examples

``` r
is_tuning(tuning_dnn_paper())
#> [1] TRUE
```
