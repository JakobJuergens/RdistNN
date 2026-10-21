# Run one replicate of CATE estimators

Draws one dataset from the named DGP (via
[`sim_cate_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_data.md)),
then runs each requested method on the **same** dataset and returns
per-replicate performance metrics. The paper's study-specific workflow
is in the installed `paper-replication/affine` bundle; these are general
simulation tools.

## Usage

``` r
sim_cate_one_rep(
  n,
  kx,
  k,
  s,
  design,
  methods = "dnn_oracle",
  tuning = NULL,
  tuning_mode = c("default", "grid"),
  alpha = 0.05,
  seed = NULL,
  sigma = 0.5,
  clip = 0.001,
  n_folds = 5L
)
```

## Arguments

- n:

  Integer sample size.

- kx:

  Integer covariate dimension (\\\ge 2\\).

- k:

  Integer DNN neighbour count.

- s:

  Integer kernel order (undersmoothing parameter). May be omitted when
  `tuning` is supplied.

- design:

  Character DGP name; passed to
  [`sim_cate_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_data.md).

- methods:

  Character vector of method names; see
  [`sim_cate_methods()`](https://JakobJuergens.github.io/RdistNN/reference/sim_cate_methods.md).
  Default: `"dnn_oracle"`.

- tuning:

  Optional kernel-order tuning specification.

- tuning_mode:

  One of `"default"` or `"grid"`. `"default"` evaluates one scalar
  kernel order. `"grid"` requires `tuning` and emits one row per DNN
  candidate while non-DNN benchmarks remain one row per replicate.

- alpha:

  Nominal interval noncoverage probability (default `0.05`, giving
  95-percent confidence intervals).

- seed:

  Optional integer seed for tuning, data and fitting, used locally
  without changing the caller's saved RNG state. `NULL` uses the current
  stream.

- sigma:

  Noise SD forwarded to the DGP.

- clip:

  Propensity clipping bound forwarded to the DGP and cross-fitting.

- n_folds:

  Number of cross-fitting folds for the `dnn_crossfit` method.

## Value

A `data.frame` with one or more rows per method and columns `method`,
`candidate_id`, `lambda`, `s`, `tau_hat`, `tau0`, `bias`, `sq_err`,
`se`, `cover` (0/1), `ci_len`, `neff`, and `status`. Methods based on
GRF also record `method_variant`, `forest_nuisance`, and whether `Y.hat`
and `W.hat` were supplied. The maintained `causal_forest` route uses
externally supplied oracle nuisance inputs. Methods with explicit
localization weights also report the realized DNN-weighted truth,
realized localization drift, estimator remainder, studentized error,
physical 90-percent-mass radius, treatment-arm effective sample sizes,
and the top-five score-contribution share. The `localized_truth` field
is the finite-sample object \\\sum_i W\_{i,s}(x_0)\tau_0(X_i)\\, not the
deterministic population branch. With \\x_0\\ and \\(n,s)\\ fixed,
expectation over the i.i.d. covariate sample and complete-DNN
exchangeability give \$\$E(\mathrm{localized\\truth})=\alpha_s(x_0).\$\$
The `localization_drift` field stores realized localization drift. The
estimator remainder combines centered score variation with fitted-state
effects. Failed methods have `NA` for all metric columns and
`status = "failed"`; registered route metadata is retained when it is
known from the selected method.

## Examples

``` r
rep1 <- sim_cate_one_rep(
  n = 200, kx = 2, k = 1, s = 12,
  design = "dgp2_nonlinear",
  methods = c("dnn_fullsample", "ols_interact"),
  seed = 1
)
rep1
#>           method candidate_id lambda  s   tau_hat     tau0       bias
#> 1 dnn_fullsample         <NA>     NA 12 1.0398647 1.248121 -0.2082565
#> 2   ols_interact         <NA>     NA NA 0.7781338 1.248121 -0.4699874
#>       sq_err        se cover    ci_len     neff localized_truth
#> 1 0.04337079 0.1810111     1 0.7095503 31.06479        1.115052
#> 2 0.22088815 0.1239821     0 0.4860010       NA              NA
#>   localization_drift estimator_remainder studentized_error score_sd0
#> 1         -0.1330695         -0.07518702         -1.150518  1.031413
#> 2                 NA                  NA                NA        NA
#>   standardized_localization_drift physical_radius_90 treated_ess control_ess
#> 1                      -0.7190849          0.7648358    20.04691    11.14627
#> 2                              NA                 NA          NA          NA
#>   top5_score_share method_variant forest_nuisance y_hat_supplied w_hat_supplied
#> 1        0.3774768           <NA>            <NA>             NA             NA
#> 2               NA           <NA>            <NA>             NA             NA
#>   status
#> 1     ok
#> 2     ok
```
