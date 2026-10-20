# Draw one sample from a named CATE simulation DGP

Generates a synthetic dataset for the CATE simulation study, including
oracle nuisance quantities \\e\\, \\m_0\\, \\m_1\\, the AIPW score used
by the general DNN CATE benchmarks.

## Usage

``` r
sim_cate_data(
  design,
  n,
  kx = 2L,
  sigma = 0.5,
  clip = 0.001,
  seed = NULL,
  x0 = NULL,
  ...
)
```

## Arguments

- design:

  Character; one of `"dgp1_linear"`, `"dgp2_nonlinear"`,
  `"dgp3_sparse_overlap"`, `"dgp4_clustered"`, `"dgp5_correlated"`, or
  `"dgp6_nuisance_efficiency"`.

- n:

  Integer sample size.

- kx:

  Integer covariate dimension (must be \\\ge 2\\).

- sigma:

  Noise standard deviation in \\Y\\.

- clip:

  Propensity clipping bound in \\(0, 0.5)\\.

- seed:

  Optional integer seed for reproducibility.

- x0:

  Optional numeric length-`kx` evaluation point. Defaults are
  DGP-specific interior points: \\(0,\ldots,0)\\ for `dgp1_linear`,
  \\(1,0,\ldots)\\ for `dgp2_nonlinear`, `dgp3_sparse_overlap`, and
  `dgp4_clustered`, and `dgp6_nuisance_efficiency`, and
  \\(0.75,0,\ldots)\\ for `dgp5_correlated`.

- ...:

  Extra arguments forwarded to the underlying DGP (e.g. `r` and
  `c_sharp` for `dgp3_sparse_overlap`).

## Value

A list with:

- `data`: a `data.frame` with columns `X1..Xkx`, `D`, `Y`, and oracle
  nuisances `e_oracle`, `m_oracle`, `m0_oracle`, `m1_oracle`,
  `tau_oracle`, and `aipw_oracle`.

- `x0`: numeric length-`kx` evaluation point.

- `tau_x0`: true treatment effect at `x0`.

- `truth`: named list of closures `tau_fn`, `e_fn`, `m0_fn`, `m_fn`.

## Details

Six designs are available:

- `dgp1_linear`:

  Linear CATE \\\tau(x) = \beta'x\\; logistic propensity. Correctly
  specified for the OLS-interaction benchmark.

- `dgp2_nonlinear`:

  Sharpened nonlinear CATE \\\tau(x) = 0.5 + 0.75\sin(1.5 x_1)\\; same
  propensity. Misspecified for OLS; primary horse-race DGP. Default
  target is the interior point \\x_0 = (1,0,\ldots)\\ so the local CATE
  is visibly nonlinear around the evaluation point.

- `dgp3_sparse_overlap`:

  Same \\\tau\\ as DGP 2; propensity interpolates between 0.5 at `x0`
  and a much sharper near-boundary logistic far from `x0`, using a
  tighter Gaussian overlap pocket. Tests estimators that do not adapt to
  local effective sample size. Default target is the same interior point
  \\x_0 = (1,0,\ldots)\\.

- `dgp4_clustered`:

  Same \\\tau\\ and propensity as DGP 2; covariates from a
  three-component mixture: 80\\ `x0`) plus two diffuse 10\\
  \\(4.5,2,0,\ldots)\\ and \\(-2.5,-2,0,\ldots)\\ with covariance
  \\1.25^2 I\\. Tests density adaptivity: the DNN s-NN bandwidth
  automatically shrinks at the dense `x0`, while a fixed-bandwidth
  kernel over-smooths. Default target is the interior point \\x_0 =
  (1,0,\ldots)\\ and the inner cluster is recentered there.

- `dgp5_correlated`:

  Covariates \\X \sim N(0,\Sigma)\\ with \\\Sigma\_{12} = \rho = 0.85\\;
  \\\tau(x) = \sin(2x_1) + 0.5x_1^2\\; evaluation point at \\x_0 =
  (0.75, 0)\\. Tests shape adaptivity: L2 and diagonal-scaled
  neighbourhoods are misaligned with the elliptical covariate
  isocontours; the `dnn_oracle_mahal` benchmark (Mahalanobis distance)
  corrects for the off-diagonal covariance.

- `dgp6_nuisance_efficiency`:

  Propensity \\e(x) = \operatorname{logit}^{-1}(1.25(x_1 - 1) -
  0.75x_2)\\ and marginal regression \\m(x) = 1 + 1.5(x_1 - 1) - x_2\\
  are aligned with simple benchmark nuisance models, while \\\tau(x) =
  0.5 + 0.60\sin(1.5 x_1)\\ remains nonlinear. Tests the precision gap
  between `dnn_fullsample` and `dnn_crossfit`. Default target is \\x_0 =
  (1,0,\ldots)\\.

## Examples

``` r
dgp <- sim_cate_data("dgp2_nonlinear", n = 200, kx = 2, seed = 1)
head(dgp$data)
#>           X1         X2 D          Y  e_oracle    m_oracle   m0_oracle
#> 1 -0.6264538  0.4094018 1 -0.1027314 0.3733369 -0.07525289 -0.03585565
#> 2  0.1836433  1.6888733 0 -0.6536643 0.3202518 -0.27977405 -0.50522996
#> 3 -0.8356286  1.5865884 1 -2.2263866 0.2295050 -1.47099798 -1.42221705
#> 4  1.5952808 -0.3309078 0  2.9318349 0.7237409  3.65752847  2.92618860
#> 5  0.3295078 -2.2852355 1  4.9663302 0.7870730  4.28830964  3.61474331
#> 6 -0.8204684  2.4976616 0 -1.5211462 0.1598876 -2.35123350 -2.31812997
#>    m1_oracle tau_oracle  aipw_oracle
#> 1 -0.1413830 -0.1055273 -0.001997335
#> 2  0.1987659  0.7039958  0.922362435
#> 3 -1.6347655 -0.2125485 -2.790362454
#> 4  3.9366882  1.0104996  0.990060973
#> 5  4.4705297  0.8557864  1.485715869
#> 6 -2.5251725 -0.2070425 -1.155705703
dgp$tau_x0
#> [1] 1.248121
```
