# Register a distance spec under a name

Register a distance spec under a name

## Usage

``` r
register_distance(name, spec)
```

## Arguments

- name:

  Character scalar registry key.

- spec:

  Either a distance spec or a function `(X, x0, ctx) -> numeric(n)`.

## Value

(Invisibly) the registered spec.

## Examples

``` r
nm <- "tmp_l1_demo"
register_distance(
  nm,
  function(X, x0, ctx = NULL) {
    rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  }
)
nm %in% distance_names()
#> [1] TRUE
```
