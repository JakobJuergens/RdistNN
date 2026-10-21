# Minimal nonparametric regression demo data (oracle m and Y only)

Generates a simple nonparametric regression dataset with features \\X\\
and response \\Y = m(X) + \epsilon\\, returning the oracle regression
function \\m(x)\\. This is designed for experimenting with the Y-only
singleton score
[`score_np_y()`](https://JakobJuergens.github.io/RdistNN/reference/score_np_y.md).

## Usage

``` r
np_demo_data(
  n = 200,
  p = 2,
  design = c("linear", "sinusoidal", "bump", "constant"),
  sigma = 0.5,
  seed = NULL,
  x0 = NULL,
  inject = NULL,
  idx_inject = 3L
)
```

## Arguments

- n:

  Integer; sample size.

- p:

  Integer; number of features in \\X\\.

- design:

  Character; one of `"linear"`, `"sinusoidal"`, `"bump"`, `"constant"`.

- sigma:

  Noise standard deviation in \\Y\\.

- seed:

  Optional integer seed, used locally without changing the caller's RNG
  state. `NULL` uses the current stream.

- x0:

  Optional target point in \\\mathbb{R}^p\\ (defaults to the zero
  vector).

- inject:

  Optional, one of `NULL`, `"na"`, `"inf"`; if provided, inject a single
  non-finite value into `Y` at `idx_inject` (useful for tests).

- idx_inject:

  Integer index at which to inject `NA`/`Inf` (default 3).

## Value

A list with:

- `data`: `data.frame` with columns `X1..Xp`, `Y`, and `m_oracle`
  (=\\m(X)\\).

- `x0`: the target point used.

- `truth`: list with closure `m_fn(x)` for evaluating \\m(x)\\ at new
  \\X\\.

- `m_x0`: the true \\m(x_0)\\.

## Examples

``` r
demo <- np_demo_data(n = 200, p = 2, design = "sinusoidal", seed = 1)
head(demo$data)
#>           X1         X2           Y   m_oracle
#> 1 -0.6264538  0.4094018  0.28598513 -0.2512354
#> 2  0.1836433  1.6888733  1.75220705  0.8043797
#> 3 -0.8356286  1.5865884 -0.04304208  0.2584566
#> 4  1.5952808 -0.3309078  0.33179122  0.5272251
#> 5  0.3295078 -2.2852355  1.25925301  1.4673640
#> 6 -0.8204684  2.4976616  1.00601700  1.1938457
demo$m_x0
#> [1] 0
```
