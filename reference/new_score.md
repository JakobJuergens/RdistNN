# Score specification (S3) for U-stat kernels

A score (kernel) is a symmetric function evaluated on s-sized subsets.
Supply `eval` or `eval_batch` for
[`score_eval()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval.md)
and
[`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md).

## Usage

``` r
new_score(
  name,
  order,
  setup = NULL,
  eval = NULL,
  eval_batch = NULL,
  eval_param = NULL,
  eval_param_batch = NULL,
  grad = NULL,
  grad_batch = NULL,
  fd_step = 1e-06,
  fd_scheme = "central",
  ptr = NULL,
  randomized = FALSE,
  symmetric = TRUE
)
```

## Arguments

- name:

  Character scalar.

- order:

  Integer order s, or function(n, data, ctx) -\> s.

- setup:

  Optional function(data, ...) -\> ctx object (list or externalptr).

- eval:

  R fallback: function(idxs, data, ctx, rng) -\> numeric(1).

- eval_batch:

  R fast path: function(subs, data, ctx, rng) -\> numeric(nrow(subs)).

- eval_param:

  Optional parameter-aware scalar evaluator: function(idxs, data, theta,
  ctx, rng) -\> numeric(1).

- eval_param_batch:

  Optional parameter-aware batch evaluator: function(subs, data, theta,
  ctx, rng) -\> numeric(nrow(subs)).

- grad:

  Optional analytic scalar gradient: function(idxs, data, theta, ctx,
  rng) -\> numeric(length(theta)).

- grad_batch:

  Optional analytic batch gradient: function(subs, data, theta, ctx,
  rng) -\> matrix(nrow(subs), length(theta)).

- fd_step:

  Finite-difference step used by
  [`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md)
  when `method = "fd"` and no explicit step is supplied.

- fd_scheme:

  Finite-difference scheme used by
  [`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md)
  when `method = "fd"` and no explicit scheme is supplied. One of
  `c("central", "forward")`.

- ptr:

  Optional compiled-score pointer metadata. The public R evaluators do
  not dispatch through this pointer; they require `eval` or
  `eval_batch`.

- randomized:

  Logical: does the kernel draw randomness?

- symmetric:

  Logical: permutation symmetry (engine assumes TRUE).

## Value

An object of class `rdistnn_ustat_score`.

## Details

If your score depends on a parameter vector (for example Z-estimation),
use `eval_param`/`eval_param_batch` and optionally `grad`/`grad_batch`
to enable derivative-aware workflows via
[`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md).

Scalar evaluators return one finite numeric value. The scalar-to-batch
fallback checks each value before assembly; logical values and multiple
values are rejected rather than coerced or truncated. Analytic batch
gradients must have exactly two dimensions: subsets by parameter
components. For a scalar parameter, a numeric vector with one entry per
subset is also accepted and converted to a one-column matrix.

## See also

[`as_score()`](https://JakobJuergens.github.io/RdistNN/reference/as_score.md),
[`score_eval()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval.md),
[`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md),
[`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md)

## Examples

``` r
sc <- new_score(
  name = "sum2",
  order = 2L,
  eval = function(idxs, data, ctx, rng) sum(data[idxs])
)
is_score(sc)
#> [1] TRUE
```
