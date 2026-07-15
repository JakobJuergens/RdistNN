# Draw one sample from a named IV-GMM simulation DGP

Generates a synthetic local IV-GMM dataset with oracle nuisance
quantities for the doubly-robust IV moment \$\$ (Z_k - \pi_k(X))\\Y -
\mu_Y(X) - (D - \mu_D(X))\theta\\. \$\$ The returned data are
estimator-free and are intended for simulation benchmarks, examples, and
future
[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
integration.

## Usage

``` r
sim_iv_data(
  design,
  n,
  kx = 2L,
  q = 2L,
  sigma = 0.5,
  clip = 0.001,
  seed = NULL,
  x0 = NULL,
  ...
)
```

## Arguments

- design:

  Character DGP name. One of `"dgp1_iv_justid_linear"`,
  `"dgp2_iv_homog_overid"`, `"dgp3_iv_heterog_overid"`,
  `"dgp4_iv_weak_local"`, `"dgp5_iv_correlated_instruments"`, or
  `"dgp6_iv_nonlinear_firststage"`.

- n:

  Integer sample size.

- kx:

  Integer covariate dimension, at least 2.

- q:

  Integer number of instruments. Phase-A simulation DGPs support
  `1 <= q <= 3`. The just-identified design always returns `q = 1`.

- sigma:

  Noise standard deviation in `Y`.

- clip:

  Clipping bound in `(0, 0.5)` for binary instrument probabilities.

- seed:

  Optional integer seed for reproducibility.

- x0:

  Optional numeric length-`kx` evaluation point. When omitted, a
  design-specific interior point is used.

- ...:

  Extra arguments forwarded to the underlying DGP. Currently `rho` is
  used by `"dgp5_iv_correlated_instruments"`.

## Value

A list with:

- `data`: a `data.frame` with `X1..Xkx`, `Y`, `D`, `Z1..Zq`, oracle
  nuisance columns `piZ*_oracle`, `muY_oracle`, `muD_oracle`, and the
  target column `theta_oracle`.

- `x0`: numeric length-`kx` evaluation point.

- `q`: realized number of instruments.

- `theta_x0`: true local IV target at `x0`.

- `truth`: named list of oracle functions.

## Examples

``` r
dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 200, kx = 2, q = 2, seed = 1)
head(dgp$data)
#>           X1         X2          Y D Z1 Z2 piZ1_oracle piZ2_oracle muY_oracle
#> 1 -0.6264538  0.4094018 -0.4726267 0  1  0   0.3842354   0.4325121  0.1585522
#> 2  0.1836433  1.6888733  2.4378363 1  0  1   0.3756758   0.4236169  0.5734972
#> 3 -0.8356286  1.5865884 -0.3643560 1  1  0   0.3017208   0.3454461 -0.6673652
#> 4  1.5952808 -0.3309078  2.5523677 0  0  0   0.6203597   0.6662061  3.5175028
#> 5  0.3295078 -2.2852355  2.8350715 1  1  1   0.6310297   0.6762599  2.9869289
#> 6 -0.8204684  2.4976616 -1.8084327 0  0  0   0.2570078   0.2970095 -1.0987363
#>   muD_oracle theta_oracle
#> 1  0.3538041   -0.0290925
#> 2  0.3493445    0.6706575
#> 3  0.2588385   -0.1485186
#> 4  0.6392690    1.1838023
#> 5  0.6433009    0.8002528
#> 6  0.2083843   -0.1412635
dgp$theta_x0
#> [1] 1.211738
```
