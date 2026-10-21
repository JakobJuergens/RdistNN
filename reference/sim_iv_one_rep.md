# Run one replicate of IV-GMM estimators

Draws one dataset from a named IV DGP via
[`sim_iv_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_data.md),
then evaluates all requested methods on that same draw.

## Usage

``` r
sim_iv_one_rep(
  n,
  kx,
  q,
  k,
  s,
  design,
  methods = "dnn_iv_oracle",
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

  Integer covariate dimension.

- q:

  Integer number of instruments requested from the DGP.

- k, s:

  Integer DNN tuning parameters.

- design:

  IV DGP name; passed to
  [`sim_iv_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_data.md).

- methods:

  Character vector of method names; see
  [`sim_iv_methods()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_methods.md).

- alpha:

  Nominal test size for confidence intervals and J-tests.

- seed:

  Optional integer seed for this replicate, used locally without
  changing the caller's saved RNG state. `NULL` uses the current stream.

- sigma:

  Noise SD forwarded to
  [`sim_iv_data()`](https://JakobJuergens.github.io/RdistNN/reference/sim_iv_data.md).

- clip:

  Clipping bound forwarded to DGPs and simulation nuisance fits.

- n_folds:

  Number of cross-fitting folds for cross-fit methods.

## Value

A `data.frame` with one row per method and columns `method`,
`theta_hat`, `theta0`, `bias`, `sq_err`, `se`, `cover`, `ci_len`,
`neff`, `J`, `df`, `pval_J`, and `reject_J`.

## Examples

``` r
sim_iv_one_rep(
  n = 120, kx = 2, q = 2, k = 1, s = 10,
  design = "dgp2_iv_homog_overid",
  methods = c("dnn_iv_oracle", "iv_oracle"),
  seed = 1
)
#>          method theta_hat   theta0        bias       sq_err        se cover
#> 1 dnn_iv_oracle  1.196561 1.211738 -0.01517726 0.0002303491 0.3300773     1
#> 2     iv_oracle  1.196561 1.211738 -0.01517726 0.0002303492 0.3300773     1
#>     ci_len     neff       J df    pval_J reject_J
#> 1 1.293879 21.94304 1.66955  1 0.1963189        0
#> 2 1.293879 21.94304 1.66955  1 0.1963188        0
```
