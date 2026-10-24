# Build a modular differentiable path for a score

Returns closures for value and gradient evaluation of a score at any
parameter vector `theta`, using analytic or finite-difference
derivatives. The `data` argument is evaluated when the path is
constructed.

## Usage

``` r
score_diff_path(
  score,
  data,
  ctx = NULL,
  method = c("auto", "analytic", "fd"),
  fd_step = NULL,
  fd_scheme = NULL
)
```

## Arguments

- score:

  A `rdistnn_ustat_score`.

- data:

  Data object passed to score evaluators.

- ctx:

  Optional context from
  [`score_setup()`](https://JakobJuergens.github.io/RdistNN/reference/score_setup.md).
  If `NULL`, one is created via
  [`score_setup()`](https://JakobJuergens.github.io/RdistNN/reference/score_setup.md).

- method:

  One of `c("auto", "analytic", "fd")` for gradients.

- fd_step:

  Optional finite-difference step override.

- fd_scheme:

  Optional finite-difference scheme override.

## Value

A list with closures:

- `value(idxs, theta, rng = NULL)` for one subset;

- `value_batch(subs, theta, rng = NULL)` for many subsets;

- `grad(idxs, theta, rng = NULL)` for one-subset gradients;

- `grad_batch(subs, theta, rng = NULL)` for batch gradients;

- `meta` with derivative settings.

## See also

[`new_score()`](https://JakobJuergens.github.io/RdistNN/reference/new_score.md),
[`score_setup()`](https://JakobJuergens.github.io/RdistNN/reference/score_setup.md),
[`score_eval()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval.md),
[`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md).

## Examples

``` r
sc <- new_score(
  name = "affine_order1",
  order = 1L,
  eval_param = function(idxs, data, theta, ctx, rng) {
    theta[1] * data[idxs] + theta[2]
  }
)
path <- score_diff_path(sc, data = c(1, 2, 4), method = "fd")
path$value(2L, theta = c(0.5, 1))
#> [1] 2
path$grad(2L, theta = c(0.5, 1))
#> [1] 2 1
```
