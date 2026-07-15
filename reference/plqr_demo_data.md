# Minimal PLQR demo data (oracle nuisances; no internal fitting)

Generates data for a partially linear quantile regression (PLQR) model
\$\$ Q\_{\tau}(Y \mid D, U) \\=\\ D^{\top}\\\theta\_{0} \\+\\ g\_{0}(U)
\\ , \$\$ with a user-chosen design for the nonparametric shape nuisance
\\g\_{0}\\; the demo specializes to \\U = X\\. Errors follow an
Asymmetric Laplace distribution (ALD) at level \\\tau\\, so the
conditional \\\tau\\-quantile equals the right-hand side.

## Usage

``` r
plqr_demo_data(
  n = 300,
  p = 2,
  p_d = 1,
  tau = 0.5,
  sigma = 0.5,
  design = c("linear", "sinusoidal", "bump"),
  seed = NULL,
  x0 = NULL,
  inject = NULL,
  inject_col = c("Y", "g_oracle", "D1", "mD1", "rD1")[1],
  idx_inject = 3L
)
```

## Arguments

- n:

  Integer; sample size.

- p:

  Integer; number of features in \\X\\ (recommend \\p \ge 2\\).

- p_d:

  Integer; dimension of \\D\\ (number of partially linear regressors).

- tau:

  Quantile level in \\(0,1)\\.

- sigma:

  Positive ALD scale for the regression error (larger means noisier
  \\Y\\).

- design:

  Character; one of `"linear"`, `"sinusoidal"`, `"bump"` for \\g\_{0}\\.

- seed:

  Optional integer seed.

- x0:

  Optional target point in \\\mathbb{R}^{p}\\ for localization (defaults
  to the zero vector).

- inject:

  Optional, one of `NULL`, `"na"`, `"inf"`; if provided, inject a single
  non-finite into `inject_col` at row `idx_inject`. Useful for testing
  runtime guards.

- inject_col:

  Character; which column to corrupt when `inject` is not `NULL`. One of
  `"Y"`, `"g_oracle"`, `"D1"`, `"mD1"`, `"rD1"` (extend as needed).

- idx_inject:

  Integer row index to corrupt (default 3).

## Value

A list with:

- `data`: a `data.frame` with columns `X1..Xp`, `D1..Dp_d`, `Y`; oracle
  nuisances `g_oracle` (the value of \\g\_{0}(X)\\), `mD1..mDp_d`
  (components of \\m\_{D}(X)\\), residuals `rD1..rDp_d` (components of
  \\D - m\_{D}(X)\\), and `density_oracle`, the true conditional density
  of \\Y\\ at the quantile threshold. Because the demo's ALD errors are
  independent of \\(X, D)\\ and homoskedastic, this density is constant,
  \\f(0) = 4\tau(1-\tau)/\sigma\\.

- `x0`: the target point used.

- `truth`: a list with `theta_true` (length `p_d`), `g_fn(x)` and
  `mD_fn(x)` closures to evaluate ground truth at new \\X\\, and
  `density_fn(v)`, the error density at signed distance \\v\\ from the
  quantile threshold (`density_fn(0)` is the `density_oracle` value).

- `g_x0`: the true value \\g\_{0}(x\_{0})\\.

## Details

The function also returns the oracle \\m\_{D}(U) \equiv \mathbb{E}\[D
\mid U\]\\ to build the orthogonal PLQR score \$\$ \psi(W;\theta) \\=\\
(\tau - \mathbf{1}\\\\Y \le D^{\top}\theta + g\_{0}(U)\\\\)\\ (D -
m\_{D}(U))\\ . \$\$

## Examples

``` r
demo <- plqr_demo_data(
  n = 300, p = 2, p_d = 1, tau = 0.5,
  design = "sinusoidal", seed = 1
)
head(demo$data)
#>           X1         X2         D1          Y    g_oracle        mD1
#> 1 -0.6264538  0.8936737 -0.9379820 -1.0859139 -0.09347465 -0.7333418
#> 2  0.1836433 -1.0472981  1.4305600  0.9953984  0.36551477  0.5291053
#> 3 -0.8356286  1.9713374 -0.9729275  0.1560512  0.60068363 -1.2899121
#> 4  1.5952808 -0.3836321  1.4359361  1.2727712  0.53664353  1.1106213
#> 5  0.3295078  1.6541453 -0.5459575  0.6313299  0.84583783 -0.4639535
#> 6 -0.8204684  1.5122127 -1.7792064 -0.6906623  0.20596416 -1.0971661
#>           rD1 density_oracle
#> 1 -0.20464019              2
#> 2  0.90145472              2
#> 3  0.31698463              2
#> 4  0.32531481              2
#> 5 -0.08200401              2
#> 6 -0.68204031              2
demo$truth$theta_true
#> [1] 0.5
demo$g_x0
#> [1] 0
```
