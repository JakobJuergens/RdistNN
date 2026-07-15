# Evaluate a tuning spec

Evaluate a tuning spec

## Usage

``` r
tuning_eval(
  tuning,
  n,
  kx,
  k,
  se = FALSE,
  jackknife = "none",
  ctx = NULL,
  mode = c("default", "grid")
)
```

## Arguments

- tuning:

  Tuning spec or registered name.

- n:

  Integer sample size.

- kx:

  Integer localization covariate dimension.

- k:

  Integer DNN neighbor count.

- se:

  Logical; whether jackknife SEs are requested.

- jackknife:

  One of `"fast"`, `"loo"`, or `"none"`.

- ctx:

  Optional context passed to tuning evaluators.

- mode:

  One of `"default"` or `"grid"`.

## Value

- `mode = "default"`: list with `s`, `meta`, and `feasibility`.

- `mode = "grid"`: list with `candidates` data.frame and `meta`.

## Examples

``` r
tuning_eval("paper_undersmooth", n = 500, kx = 2, k = 1)
#> $s
#> [1] 42
#> 
#> $meta
#> $meta$lambda
#> [1] 0.4
#> 
#> $meta$formula
#> [1] "lambda * n^((kx + 1)/(kx + 2))"
#> 
#> $meta$rounding
#> [1] "floor"
#> 
#> $meta$name
#> [1] "paper_undersmooth"
#> 
#> $meta$n
#> [1] 500
#> 
#> $meta$kx
#> [1] 2
#> 
#> 
#> $feasibility
#> $feasibility$k_le_s
#> [1] TRUE
#> 
#> $feasibility$s_lt_n
#> [1] TRUE
#> 
#> 
tuning_eval("paper_undersmooth", n = 500, kx = 2, k = 1, mode = "grid")
#> $candidates
#>   lambda  s feasible reason
#> 1   0.30 31     TRUE     ok
#> 2   0.40 42     TRUE     ok
#> 3   0.50 52     TRUE     ok
#> 4   0.75 79     TRUE     ok
#> 
#> $meta
#> $meta$name
#> [1] "paper_undersmooth"
#> 
#> $meta$n
#> [1] 500
#> 
#> $meta$kx
#> [1] 2
#> 
#> 
```
