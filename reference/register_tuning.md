# Register a tuning spec under a name

Register a tuning spec under a name

## Usage

``` r
register_tuning(name, tuning)
```

## Arguments

- name:

  Character scalar registry key.

- tuning:

  A tuning spec or function coercible via
  [`as_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/as_tuning.md).

## Value

(Invisibly) the registered spec.

## Examples

``` r
register_tuning("tmp_fixed_4", function(...) 4L)
"tmp_fixed_4" %in% tuning_names()
#> [1] TRUE
```
