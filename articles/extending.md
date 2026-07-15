# Extending RdistNN: custom distances, rank weights, and scores

## Why The Package Is Modular

`RdistNN` is built in layers:

- ranking: turn `X` and `x0` / `X0` into distances, orders, and ranks;
- weighting: turn ranks into localization weights;
- tuning: resolve kernel order `s` from `(n, kx, k)` using a pluggable
  rule;
- scoring: evaluate a kernel or estimating function on subsets of the
  data;
- wrappers: combine the pieces into user-facing workflows such as
  [`np_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/np_dnn.md),
  [`cate_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/cate_dnn.md),
  or
  [`plqr_dnn()`](https://JakobJuergens.github.io/RdistNN/reference/plqr_dnn.md).

That separation is deliberate. Most new methods only need to replace one
layer.

As a rule of thumb:

- use
  [`as_distance()`](https://JakobJuergens.github.io/RdistNN/reference/as_distance.md),
  [`as_rank_weight()`](https://JakobJuergens.github.io/RdistNN/reference/as_rank_weight.md),
  or
  [`as_score()`](https://JakobJuergens.github.io/RdistNN/reference/as_score.md)
  for quick scalar prototypes;
- use
  [`new_distance()`](https://JakobJuergens.github.io/RdistNN/reference/new_distance.md),
  [`new_rank_weight()`](https://JakobJuergens.github.io/RdistNN/reference/new_rank_weight.md),
  or
  [`new_score()`](https://JakobJuergens.github.io/RdistNN/reference/new_score.md)
  when you need explicit batch evaluation, reusable setup context,
  registration, or parameter-aware score derivatives.

## Choose The Right Extension Point

- If you want a different notion of closeness, implement a custom
  distance or a custom ranking provider.
- If you like the ranks but want a different localization profile,
  implement a custom rank-weight rule.
- If you want a different kernel-order rule, implement a custom tuning
  spec.
- If the localization stays fixed but the target moment changes,
  implement a custom score.
- If you already have distances, orders, or ranks from another engine,
  pass them directly into
  [`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md)
  instead of rebuilding that layer inside `RdistNN`.

## 0) Kernel-Order Tuning Specs

Tuning specs mirror the same constructor/registry style as distances and
rank weights. The built-in paper default is:

``` r

tun <- tuning_dnn_paper(lambda = 1)
tuning_eval(tun, n = 500, kx = 2, k = 1)
tuning_eval(tun, n = 500, kx = 2, k = 1, mode = "grid")
```

You can add custom rules:

``` r

tun_fast <- as_tuning(function(n, kx, ...) floor(0.8 * n ^ ((kx + 1) / (kx + 2))), name = "fast_rule")
register_tuning("fast_rule", tun_fast)
tuning_names()
```

## 1) Custom Distances

### Quick Prototype With `as_distance()`

Use
[`as_distance()`](https://JakobJuergens.github.io/RdistNN/reference/as_distance.md)
when a simple function `(X, x0, ctx) -> numeric(n)` is enough.

``` r

# L-infinity distance
linf <- as_distance(
  function(X, x0, ctx = NULL) {
    apply(abs(X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)), 1L, max)
  },
  name = "linf"
)

register_distance("linf", linf)
distance_names()
```

You can test the distance directly before registering it:

``` r

X <- matrix(rnorm(40), 20, 2)
x0 <- c(0, 0)

distance_eval(linf, X, x0)
```

### Reusable Distance With `new_distance()`, `setup()`, and `eval_batch()`

Use
[`new_distance()`](https://JakobJuergens.github.io/RdistNN/reference/new_distance.md)
when you want to precompute context once on `X` and reuse it across many
targets.

``` r

scaled_l2 <- new_distance(
  name = "scaled_l2_demo",
  setup = function(X) {
    scale <- apply(X, 2L, stats::sd)
    scale[!is.finite(scale) | scale <= 0] <- 1
    list(scale = scale)
  },
  eval_batch = function(X, X0, ctx) {
    Xs <- sweep(X, 2L, ctx$scale, "/")
    X0s <- sweep(X0, 2L, ctx$scale, "/")

    vapply(seq_len(nrow(X0s)), function(j) {
      rowSums((Xs - matrix(X0s[j, ], nrow(Xs), ncol(Xs), byrow = TRUE))^2)
    }, numeric(nrow(Xs)))
  }
)

X <- matrix(rnorm(80), 40, 2)
X0 <- matrix(rnorm(12), 6, 2)

ctx_dist <- scaled_l2$setup(X)
D <- distance_eval_batch(scaled_l2, X, X0, ctx = ctx_dist)
dim(D)
```

When you call
[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md)
with a distance spec and `ctx = NULL`, `RdistNN` automatically uses the
distance-side `setup(X)` if it exists.

## 2) Custom Rank-Weight Rules

Rank weights only see integer ranks, not raw distances. Rank `1` means
“closest observation,” rank `2` means “second closest,” and so on.

That separation is useful because you can change the localization
profile without changing the ranking logic.

### Quick Prototype With `as_rank_weight()`

``` r

rw_inv <- as_rank_weight(
  function(ranks, ctx = NULL) {
    w <- 1 / as.numeric(ranks)
    w / sum(w)
  },
  name = "inv_rank"
)

register_rank_weight("inv_rank", rw_inv)
rank_weight_names()
```

### A Compact-Support Rule

This example keeps only the first `m` closest observations in rank
order. That is a generic finite-support rule. It is not the same thing
as exact DNN.

``` r

make_top_m <- function(m) {
  stopifnot(length(m) == 1L, m >= 1L)

  as_rank_weight(
    function(ranks) {
      keep <- as.numeric(ranks <= m)
      keep / sum(keep)
    },
    name = paste0("top_", m)
  )
}

rw_top10 <- make_top_m(10L)
```

### Batched Rank Weights and Shared `ctx`

For performance, prefer `eval_batch()` once the prototype is stable.

``` r

rw_exp <- new_rank_weight(
  name = "exp_decay",
  eval_batch = function(R, ctx = NULL) {
    lambda <- if (is.null(ctx) || is.null(ctx$lambda)) 0.25 else ctx$lambda
    W <- exp(-lambda * (R - 1L))
    sweep(W, 2L, colSums(W), "/")
  }
)

R <- cbind(1:5, c(2L, 1L, 3L, 4L, 5L))
rank_weight_eval_batch(rw_exp, R, ctx = list(lambda = 0.5))
```

At the top level,
[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md)
passes any supplied `ctx` through to both the ranking and weight
evaluators. If both sides need context, build one shared list.

Also note one current API detail: `new_rank_weight(..., setup = ...)` is
a low-level hook, but
[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md)
does not automatically call rank-weight-side `setup()` for you. In
practice, the cleanest high-level patterns are:

- close over fixed tuning parameters inside the constructor; or
- build a shared `ctx` object yourself and pass it into
  [`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md).

## 3) Put Ranking And Weighting Together With `weights_matrix()`

[`weights_matrix()`](https://JakobJuergens.github.io/RdistNN/reference/weights_matrix.md)
is the main modular entry point. The `ranking` argument can be:

- a distance spec or registered distance name;
- a function `(X, X0, ctx)` that returns `dx`, `ord`, or `ranks`;
- a precomputed list with `dx`, `ord`, or `ranks`;
- a ranks matrix directly.

That makes it the right place to integrate custom ranking code,
precomputed distances, or external nearest-neighbor engines.

``` r

n <- 80
p <- 2
X <- matrix(rnorm(n * p), n, p)
X0 <- matrix(rnorm(4 * p), 4, p)

# If only the distance needs setup, leave ctx = NULL and let the distance
# auto-setup path run.
out1 <- weights_matrix(
  X = X,
  X0 = X0,
  ranking = scaled_l2,
  weight = rw_top10,
  return = c("W", "dx", "ord", "ranks")
)

# If both the distance and the weight need context, build one shared list.
ctx_both <- c(scaled_l2$setup(X), list(lambda = 0.35))
out2 <- weights_matrix(
  X = X,
  X0 = X0,
  ranking = scaled_l2,
  weight = rw_exp,
  ctx = ctx_both,
  return = c("W", "dx", "ord", "ranks")
)

dim(out2$W)
```

For debugging custom ranking behavior, request the intermediate objects:

``` r

out_dbg <- weights_matrix(
  X = X,
  X0 = X0,
  ranking = scaled_l2,
  weight = rw_exp,
  ctx = ctx_both,
  return = c("dx", "ord", "ranks", "W")
)

head(out_dbg$dx[, 1])
head(out_dbg$ord[, 1])
head(out_dbg$ranks[, 1])
head(out_dbg$W[, 1])
```

## 4) Custom Scores

The score layer is separate from localization. A score describes what
happens once a subset of observations has been selected.

### Quick Prototype With `as_score()`

Use
[`as_score()`](https://JakobJuergens.github.io/RdistNN/reference/as_score.md)
for a simple kernel with fixed order `s`.

``` r

sc_avg2 <- as_score(
  function(idxs, data, ctx = NULL, rng = NULL) {
    mean(data$Y[idxs])
  },
  s = 2L,
  name = "avg2"
)

register_score("avg2", sc_avg2)
list_scores()

df <- data.frame(Y = rnorm(10))
subs <- rbind(c(1L, 2L), c(4L, 8L), c(3L, 9L))
score_eval_batch(sc_avg2, subs = subs, data = df)
```

### Reusable Score With `new_score()` and `setup()`

Use `setup()` when the score repeatedly accesses the same data columns
or transformed objects.

``` r

sc_pair_gap <- new_score(
  name = "pair_gap_sq",
  order = 2L,
  setup = function(data, ...) list(y = as.numeric(data$Y)),
  eval_batch = function(subs, data, ctx, rng = NULL) {
    (ctx$y[subs[, 1L]] - ctx$y[subs[, 2L]])^2
  }
)

df <- data.frame(Y = rnorm(10))
ctx_score <- score_setup(sc_pair_gap, df)
subs <- rbind(c(1L, 2L), c(2L, 5L), c(3L, 4L))

score_eval_batch(sc_pair_gap, subs, data = df, ctx = ctx_score)
```

### Parameter-Aware Scores And Derivatives With `score_diff_path()`

If your score depends on a parameter vector, implement `eval_param` or
an `eval` that accepts `theta`, then use
[`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md)
to get reusable value and gradient closures.

``` r

sc_affine <- new_score(
  name = "affine_order1",
  order = 1L,
  eval_param = function(idxs, data, theta, ctx, rng = NULL) {
    theta[1] * data$Y[idxs] + theta[2]
  },
  grad = function(idxs, data, theta, ctx, rng = NULL) {
    c(data$Y[idxs], 1)
  }
)

df <- data.frame(Y = c(1, 2, 4))
path <- score_diff_path(sc_affine, data = df, method = "auto")

path$value(2L, theta = c(0.5, 1))
path$grad(2L, theta = c(0.5, 1))
```

If you do not provide an analytic gradient,
`score_diff_path(..., method = "fd")` uses finite differences based on
`eval_param` / `eval_param_batch`.

## 5) A Minimal Custom Estimator Pattern

The modular API is meant to let you assemble your own localized
estimators without rewriting the whole package. A minimal pattern is:

``` r

my_local_mean <- function(X, x0, y, ranking = "l2", weight) {
  W <- weights_matrix(
    X = X,
    x0 = x0,
    ranking = ranking,
    weight = weight,
    return = "W"
  )$W[, 1L]

  sum(W * y)
}
```

For more complicated estimators, replace the scalar `y` contribution
with your own score or estimating equation.

## 6) Validation And Debugging Checklist

- Start by testing the low-level evaluator directly:
  [`distance_eval()`](https://JakobJuergens.github.io/RdistNN/reference/distance_eval.md),
  [`distance_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/distance_eval_batch.md),
  [`rank_weight_eval()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_eval.md),
  [`rank_weight_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/rank_weight_eval_batch.md),
  [`score_eval()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval.md),
  or
  [`score_eval_batch()`](https://JakobJuergens.github.io/RdistNN/reference/score_eval_batch.md).
- Remember that rank weights receive ranks, not distances. Smaller rank
  means closer observation.
- Prototype with `as_*()` first. Once the logic is stable, move to
  `new_*()` with `eval_batch()` for speed.
- Use `weights_matrix(return = c("dx", "ord", "ranks", "W"))` to debug
  ranking and tie handling before worrying about the downstream
  estimator.
- If both ranking and weighting need context, create it explicitly and
  pass it through `ctx`.
- For scores, check
  [`score_order()`](https://JakobJuergens.github.io/RdistNN/reference/score_order.md)
  and
  [`score_setup()`](https://JakobJuergens.github.io/RdistNN/reference/score_setup.md)
  before wiring the score into root finding, optimization, or jackknife
  code.
- For parameter-aware scores, test
  [`score_diff_path()`](https://JakobJuergens.github.io/RdistNN/reference/score_diff_path.md)
  in isolation before embedding it in a solver.
- Register names only after the underlying object works directly. The
  registry is a convenience layer, not a debugging tool.
