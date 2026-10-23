# Kernel Order Tuning

The order `s` controls localization: larger values concentrate weight
closer to the reporting profile. A grid compares estimates across orders
and retains its specified default; it does not select an optimal order.
The
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
gives the statistical motivation and inference conditions.

## Fit and compare a grid

``` r

demo <- np_demo_data(n = 200, p = 2, seed = 1)
fit <- np_dnn(
  data = demo$data, x_cols = c("X1", "X2"), y_col = "Y",
  x0 = demo$x0, tuning_mode = "grid"
)
fit$tuning_grid$candidates
#>   candidate_id computed_key is_default feasible reason status lambda  s
#> 1   lambda_0p3         s_15      FALSE     TRUE     ok     ok   0.30 15
#> 2   lambda_0p4         s_21       TRUE     TRUE     ok     ok   0.40 21
#> 3   lambda_0p5         s_26      FALSE     TRUE     ok     ok   0.50 26
#> 4  lambda_0p75         s_39      FALSE     TRUE     ok     ok   0.75 39
fit$tuning_grid$estimates
#>   candidate_id computed_key lambda  s target_id component  estimate         se
#> 1   lambda_0p3         s_15   0.30 15         1      <NA> 0.1311178 0.08214191
#> 2   lambda_0p4         s_21   0.40 21         1      <NA> 0.1650251 0.08880114
#> 3   lambda_0p5         s_26   0.50 26         1      <NA> 0.1835091 0.09559482
#> 4  lambda_0p75         s_39   0.75 39         1      <NA> 0.2118122 0.11896597
#>   status reason
#> 1     ok     ok
#> 2     ok     ok
#> 3     ok     ok
#> 4     ok     ok
```

The main estimate remains the default candidate. Use a candidate view to
inspect another order without fitting the nuisance functions again:

``` r

view <- candidate_fit(fit, 1L)
coef(view)
#>  target_1 
#> 0.1311178
```

Candidate identifiers refer to rows of the tuning grid, not just the
resolved order. If two tuning settings resolve to the same `s`, they
remain distinct candidates even when their computation is shared. Plot
labels include their identifiers when needed to distinguish them.
Candidate reconstruction validates the retained controls, target
identities, cached ranks, per-target distance groups and selected
policy; refit an incomplete or inconsistent saved object.

CATE and CLATE use the same `tuning_mode` convention. See [Getting
started](https://JakobJuergens.github.io/RdistNN/articles/getting-started.md)
and [Conditional
LATE](https://JakobJuergens.github.io/RdistNN/articles/conditional-late.md)
for their signal inputs. `plot_tuning(fit)` draws the sensitivity
comparison when ggplot2 is installed.

## Inspect the default rule

The built-in rule is `floor(lambda * n^((kx + 1)/(kx + 2)))`, with
`lambda = 0.4` and sensitivity multipliers `c(0.3, 0.4, 0.5, 0.75)`.
Here `kx` is the number of reporting covariates. `lambda = 1` remains an
available, more concentrated localization setting. Delete-one
uncertainty requires `s < n`.

``` r

tun <- tuning_dnn_paper()
tuning_eval(tun, n = 1000, kx = 2, k = 1)$s
#> [1] 71
tuning_eval(tun, n = 1000, kx = 2, k = 1, mode = "grid")$candidates
#>   lambda   s feasible reason
#> 1   0.30  53     TRUE     ok
#> 2   0.40  71     TRUE     ok
#> 3   0.50  88     TRUE     ok
#> 4   0.75 133     TRUE     ok
```

Explicit `s` and a tuning specification are alternative inputs; grid
mode cannot be combined with explicit `s`. Diagnostics record the
resolved order and its provenance:

``` r

as.data.frame(diagnostics(fit))[, c("target_id", "tuning_name", "tuning_default_s")]
#>          target_id       tuning_name tuning_default_s
#> target_1  target_1 paper_undersmooth               21
```

## Define and register a custom rule

For a fixed order, a small function is enough:

``` r

fixed12 <- as_tuning(function(...) 12L, name = "fixed_12")
tuning_eval(fixed12, n = 300, kx = 2, k = 1)$s
#> [1] 12
register_tuning("fixed_12", fixed12)
get_tuning("fixed_12")$name
#> [1] "fixed_12"
```

[`new_tuning()`](https://JakobJuergens.github.io/RdistNN/reference/new_tuning.md)
also accepts explicit candidate generation and metadata; see its help
for the full interface. A custom rule must respect the estimator’s order
constraints. [DNN
weights](https://JakobJuergens.github.io/RdistNN/articles/dnn-weights.md)
explains the resulting rank weights, and [Notation and
conventions](https://JakobJuergens.github.io/RdistNN/articles/notation-and-conventions.md)
lists the argument and result shapes.
