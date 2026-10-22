# Create a local GMM model object

Defines a vector-moment model for use with
[`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md).
A model supplies observation-level moments and, optionally,
observation-level Jacobians and target-specific starting values.

## Usage

``` r
new_gmm_model(
  name,
  setup,
  moment_obs,
  jacobian_obs = NULL,
  start = NULL,
  coef_names,
  moment_names = NULL
)
```

## Arguments

- name:

  Non-empty model name.

- setup:

  Function that validates and prepares model data.

- moment_obs:

  Function returning an `n x q` numeric moment matrix.

- jacobian_obs:

  Optional function returning observation-level Jacobians. For scalar
  `theta`, an `n x q` matrix is accepted; otherwise use an `n x q x p`
  array.

- start:

  Optional function `(sample, ctx, w, x0)` returning a finite numeric
  starting vector, evaluated independently for each target. An explicit
  `start` passed to
  [`gmm_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/gmm_dnn.md)
  takes precedence.

- coef_names:

  Non-empty character vector naming the parameter entries.

- moment_names:

  Optional character vector naming the moment components.

## Value

An object of class `gmm_model`.

## Details

The model functions are called with named arguments:
`setup(sample, n, ...)`, `moment_obs(sample, theta, ctx)`,
`jacobian_obs(sample, theta, ctx)`, and `start(sample, ctx, w, x0)`. The
`setup()` return value is passed back as `ctx`.

## Examples

``` r
mod <- new_gmm_model(
  name = "mean_zero",
  setup = function(sample, n, ...) list(g = sample$g),
  moment_obs = function(sample, theta, ctx) matrix(ctx$g - theta[1], ncol = 1),
  jacobian_obs = function(sample, theta, ctx) matrix(-1, length(ctx$g), 1),
  start = function(sample, ctx, w, x0) sum(w * ctx$g),
  coef_names = "theta",
  moment_names = "mean"
)
is_gmm_model(mod)
#> [1] TRUE
```
