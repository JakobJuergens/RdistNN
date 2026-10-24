# Construct a tuning spec

Construct a tuning spec

## Usage

``` r
new_tuning(name, eval = NULL, eval_grid = NULL, setup = NULL)
```

## Arguments

- name:

  Character scalar.

- eval:

  Optional function `(n, kx, k, se, jackknife, ctx) -> list|numeric`.

- eval_grid:

  Optional function `(n, kx, k, se, jackknife, ctx) -> data.frame`.

- setup:

  Optional setup function for future use.

## Value

A validated `rdistnn_tuning` spec.

## Examples

``` r
tun <- new_tuning(
  name = "fixed_s",
  eval = function(n, kx, k, se, jackknife, ctx = NULL) list(s = 5L)
)
is_tuning(tun)
#> [1] TRUE
```
