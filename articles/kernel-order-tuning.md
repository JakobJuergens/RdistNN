# Kernel Order Tuning

## Overview

`RdistNN` exposes kernel-order selection through a modular tuning API:

- build specs with
  [`new_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/new_tuning.md)
  or
  [`as_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/as_tuning.md),
- register/reuse with
  [`register_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/register_tuning.md)
  and
  [`get_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/get_tuning.md),
- evaluate defaults and candidate grids via
  [`tuning_eval()`](https://JakobJuergens.github.io/RdistNN/reference/tuning_eval.md);
- inspect resolved tuning metadata in fitted-object diagnostics.

## Paper Default Undersmoothing Rule

The built-in rule is:

``` math
s(\lambda) = \left\lfloor \lambda\, n^{(k_x+1)/(k_x+2)} \right\rfloor
```

with finite-sample default `lambda = 0.4` and default sensitivity grid
`lambda_grid = c(0.3, 0.4, 0.5, 0.75)`. The value `lambda = 1` remains
available as a more aggressive localization benchmark.

``` r

tun <- tuning_dnn_paper()
tuning_eval(tun, n = 1000, kx = 2, k = 1, mode = "default")
#> $s
#> [1] 71
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
#> [1] 1000
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
```

``` r

tuning_eval(tun, n = 1000, kx = 2, k = 1, mode = "grid")$candidates
#>   lambda   s feasible reason
#> 1   0.30  53     TRUE     ok
#> 2   0.40  71     TRUE     ok
#> 3   0.50  88     TRUE     ok
#> 4   0.75 133     TRUE     ok
```

## Registry Usage

``` r

tuning_names()
#> [1] "paper_undersmooth"
get_tuning("paper_undersmooth")$name
#> [1] "paper_undersmooth"
```

## Custom Rule

``` r

fixed12 <- as_tuning(function(...) 12L, name = "fixed_12")
tuning_eval(fixed12, n = 300, kx = 2, k = 1)
#> $s
#> [1] 12
#> 
#> $meta
#> $meta$name
#> [1] "fixed_12"
#> 
#> $meta$n
#> [1] 300
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
```

A custom rule can also return richer metadata:

``` r

meta_rule <- new_tuning(
  name = "meta_rule",
  eval = function(n, kx, k, se, jackknife, ctx = NULL) {
    list(
      s = floor(0.9 * n ^ ((kx + 1) / (kx + 2))),
      meta = list(
        lambda = 0.9,
        formula = "0.9 * n^((kx + 1)/(kx + 2))",
        rounding = "floor"
      )
    )
  }
)

tuning_eval(meta_rule, n = 500, kx = 2, k = 1)
#> $s
#> [1] 95
#> 
#> $meta
#> $meta$lambda
#> [1] 0.9
#> 
#> $meta$formula
#> [1] "0.9 * n^((kx + 1)/(kx + 2))"
#> 
#> $meta$rounding
#> [1] "floor"
#> 
#> $meta$name
#> [1] "meta_rule"
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
```

## Wrapper Integration

High-level wrappers accept `tuning = ...` and resolve `s` internally:

``` r

demo <- np_demo_data(n = 200, p = 2, seed = 1)
fit <- np_dnn(
  data = demo$data,
  x_cols = c("X1", "X2"),
  y_col = "Y",
  x0 = demo$x0,
  k = 1,
  tuning = tuning_dnn_paper(),
  se = TRUE,
  return = c("W", "ess", "tuning")
)

fit$meta$s
#> [1] 21
fit$tuning$grid$candidates
#>   candidate_id computed_key is_default feasible reason status lambda  s
#> 1   lambda_0p3         s_15      FALSE     TRUE     ok     ok   0.30 15
#> 2   lambda_0p4         s_21       TRUE     TRUE     ok     ok   0.40 21
#> 3   lambda_0p5         s_26      FALSE     TRUE     ok     ok   0.50 26
#> 4  lambda_0p75         s_39      FALSE     TRUE     ok     ok   0.75 39
diagnostics(fit)[, c("target_id", "tuning_name", "tuning_default_s", "tuning_grid_n_candidates")]
#> <dnn_diagnostics>
#>   targets: 1
#>   dist   : available for 0/1 target(s)
```

## Related Articles

- [Notation and package
  conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md):
  [`vignette("notation-and-conventions", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md)
  relates `s`, `k`, `lambda`, support, and effective sample size to the
  paper notation.
- [DNN
  weights](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md):
  [`vignette("dnn-weights", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md)
  explains how the selected order enters exact and asymptotic weights.
- [Diagnostics and
  reporting](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md):
  [`vignette("diagnostics-reporting", package = "RdistNN")`](https://JakobJuergens.github.io/RdistNN/articles/diagnostics-reporting.md)
  shows how resolved tuning metadata appears in fitted-object
  diagnostics.
