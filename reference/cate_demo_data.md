# Minimal CATE demo data (oracle nuisances + pseudo-outcome)

Generates a small synthetic dataset for experimenting with the CATE
score. It returns \\X\\, \\D\\, \\Y\\, along with oracle nuisance
quantities: \$\$ m(x) \\=\\ \mathbb{E}\[Y \mid X=x\] \\=\\ m\_{0}(x)
\\+\\ e(x)\\\tau(x) \\ , \qquad e(x) \\=\\ \Pr(D=1 \mid X=x) \\ , \$\$
where \\m\_{0}(x) = \mathbb{E}\[Y \mid X=x, D=0\]\\ and \\\tau(x)\\ is
the treatment effect. It also returns the oracle AIPW pseudo-outcome.

## Usage

``` r
cate_demo_data(
  n = 200,
  p = 2,
  design = c("simple", "hetero", "bump"),
  sigma = 0.5,
  clip = 0.001,
  seed = NULL,
  x0 = NULL,
  inject = NULL,
  inject_col = "aipw_oracle",
  idx_inject = 3L
)
```

## Arguments

- n:

  Integer sample size.

- p:

  Integer feature dimension, at least 2.

- design:

  Character; one of `"simple"`, `"hetero"`, `"bump"` for \\\tau(x)\\.

- sigma:

  Noise standard deviation in \\Y\\.

- clip:

  Propensity clipping in \\(0,\\0.5)\\; applied to \\e(x)\\.

- seed:

  Optional integer seed, used locally without changing the caller's RNG
  state. `NULL` uses the current stream.

- x0:

  Optional numeric length-\\p\\ target point (defaults to the zero
  vector).

- inject:

  Optional, one of `NULL`, `"na"`, `"inf"`; if provided, inject a single
  non-finite value into one of the returned columns (see `inject_col`)
  at row `idx_inject`; useful for testing runtime guards.

- inject_col:

  Character; which column to corrupt when `inject` is not `NULL`. One of
  `"Y"`, `"aipw_oracle"`, `"e_oracle"`, `"m_oracle"`, `"m0_oracle"`,
  `"tau_oracle"`. Default `"aipw_oracle"`.

- idx_inject:

  Integer index to corrupt (default 3).

## Value

A list with:

- `data`: a `data.frame` with columns `X1..Xp`, `D`, `Y`, and oracle
  nuisances `e_oracle` (=\\e(X)\\), `m_oracle` (=\\m(X)\\), `m0_oracle`
  (=\\m\_{0}(X)\\), `m1_oracle` (treated outcome mean), `tau_oracle`
  (=\\\tau(X)\\), `aipw_oracle` (oracle AIPW signal).

- `x0`: the target point used (length \\p\\).

- `truth`: list of closures `tau_fn(x)`, `e_fn(x)`, `m0_fn(x)`,
  `m_fn(x)` to evaluate ground truth at any \\x\\.

- `tau_x0`: the true \\\tau(x\_{0})\\.

## Details

No model fitting is performed; all quantities come from the DGP.

## Examples

``` r
demo <- cate_demo_data(n = 200, p = 2, design = "hetero", seed = 1)
str(demo$data)
#> 'data.frame':    200 obs. of  10 variables:
#>  $ X1         : num  -0.626 0.184 -0.836 1.595 0.33 ...
#>  $ X2         : num  0.409 1.689 1.587 -0.331 -2.285 ...
#>  $ D          : int  1 0 1 0 1 0 0 0 1 0 ...
#>  $ Y          : num  0.21 -0.654 -1.885 2.932 4.772 ...
#>  $ e_oracle   : num  0.373 0.32 0.23 0.724 0.787 ...
#>  $ m_oracle   : num  0.0414 -0.3159 -1.3926 3.6498 4.1356 ...
#>  $ m0_oracle  : num  -0.0359 -0.5052 -1.4222 2.9262 3.6147 ...
#>  $ m1_oracle  : num  0.171 0.0861 -1.2931 3.926 4.2765 ...
#>  $ tau_oracle : num  0.207 0.591 0.129 1 0.662 ...
#>  $ aipw_oracle: num  0.31 0.81 -2.449 0.979 1.292 ...
demo$tau_x0
#> [1] 0.5
```
