# Retrieve a registered distance spec

Retrieve a registered distance spec

## Usage

``` r
get_distance(name)
```

## Arguments

- name:

  Character scalar.

## Value

A distance spec.

## Examples

``` r
d <- get_distance("l2")
is_distance(d)
#> [1] TRUE
```
