# -------------------------------------------------------------------
# ustat_score (S3) — scaffolding tests for score.R
# Tags:
#   [smoke]    Construction & printing
#   [validate] Validation errors & warnings
#   [setup]    Context plumbing
#   [order]    Order resolution (scalar vs function)
#   [engine]   Eval / Eval-batch shapes, finite checks, input-path errors
#   [helpers]  as_score / predicates
# -------------------------------------------------------------------

test_that("[smoke] can construct a minimal R score and print it", {
  f <- function(idxs, data, ctx = NULL, rng = NULL) 0.0
  sc <- new_score(name = "zero", order = 2L, eval = f)
  expect_s3_class(sc, "ustat_score")
  expect_output(print(sc), "<ustat_score:") # prefix
  expect_output(print(sc), "order: 2") # scalar order
  expect_output(print(sc), "impl :\\s*R eval") # implementation line
})

test_that("[validate] errors on empty name, bad order, missing impl; warns on non-symmetric", {
  f <- function(idxs, data, ctx = NULL, rng = NULL) 0.0

  expect_error(
    new_score(name = "", order = 2L, eval = f),
    "score name must be non-empty"
  )

  # bad order: not scalar int and not function(n, data, ctx)
  expect_error(
    new_score(name = "bad", order = "two", eval = f),
    "order must be integer or function\\(n, data, ctx\\)"
  )

  # missing implementations
  expect_error(
    new_score(name = "noimpl", order = 2L),
    "Provide at least one of: ptr, eval_batch, eval"
  )

  # non-symmetric warning
  expect_warning(
    new_score(name = "ns", order = 2L, eval = f, symmetric = FALSE),
    "non-symmetric"
  )
})

test_that("[setup] distance-style setup passes ctx into eval and eval_batch", {
  used_ctx_eval <- FALSE
  used_ctx_batch <- FALSE

  sc <- new_score(
    name = "ctx",
    order = 2L,
    setup = function(data, scale = 1) list(scale = scale),
    eval = function(idxs, data, ctx = NULL, rng = NULL) {
      used_ctx_eval <<- isTRUE(ctx$scale == 3)
      sum(data[idxs, 1]) * ctx$scale
    },
    eval_batch = function(subs, data, ctx = NULL, rng = NULL) {
      used_ctx_batch <<- isTRUE(ctx$scale == 3)
      apply(subs, 1L, function(idxs) sum(data[idxs, 1]) * ctx$scale)
    }
  )

  set.seed(1)
  data <- cbind(rnorm(5), rnorm(5))
  ctx <- score_setup(sc, data, scale = 3)

  # eval
  out1 <- score_eval(sc, c(1L, 2L), data, ctx)
  expect_true(used_ctx_eval)
  expect_equal(length(out1), 1L)
  expect_true(is.finite(out1))

  # eval_batch
  subs <- rbind(c(1L, 2L), c(3L, 4L))
  out2 <- score_eval_batch(sc, subs, data, ctx)
  expect_true(used_ctx_batch)
  expect_equal(length(out2), nrow(subs))
  expect_true(all(is.finite(out2)))
})

test_that("[order] resolves scalar order and function order correctly", {
  f <- function(idxs, data, ctx = NULL, rng = NULL) 0.0

  sc_scalar <- new_score(name = "s2", order = 2L, eval = f)
  expect_identical(score_order(sc_scalar, n = 10L, data = NULL), 2L)

  sc_fun <- new_score(
    name = "s_fun",
    order = function(n, data, ctx) as.integer(max(1, floor(n / 5))),
    eval = f
  )
  expect_identical(score_order(sc_fun, n = 10L, data = NULL), 2L)
  expect_identical(score_order(sc_fun, n = 7L, data = NULL), 1L)
})

test_that("[engine] score_eval enforces scalar finite output; score_eval_batch enforces length", {
  data <- cbind(1:5, 6:10)
  subs <- rbind(c(1L, 2L), c(3L, 4L))

  # Good eval
  sc_ok <- new_score(
    name = "ok",
    order = 2L,
    eval = function(idxs, data, ctx = NULL, rng = NULL) sum(data[idxs, 1])
  )
  o <- score_eval(sc_ok, c(1L, 2L), data)
  expect_equal(length(o), 1L)
  expect_true(is.finite(o))

  # Bad eval: vector length > 1
  sc_bad_len <- new_score(
    name = "badlen",
    order = 2L,
    eval = function(idxs, data, ctx = NULL, rng = NULL) c(1, 2)
  )
  expect_error(
    score_eval(sc_bad_len, c(1L, 2L), data),
    "must return a finite scalar numeric"
  )

  # Bad eval: non-finite
  sc_bad_nan <- new_score(
    name = "badnan",
    order = 2L,
    eval = function(idxs, data, ctx = NULL, rng = NULL) NaN
  )
  expect_error(
    score_eval(sc_bad_nan, c(1L, 2L), data),
    "must return a finite scalar numeric"
  )

  # Good eval_batch
  sc_b_ok <- new_score(
    name = "bok",
    order = 2L,
    eval_batch = function(subs, data, ctx = NULL, rng = NULL) {
      apply(subs, 1L, function(idxs) sum(data[idxs, 2]))
    }
  )
  ob <- score_eval_batch(sc_b_ok, subs, data)
  expect_length(ob, nrow(subs))
  expect_true(all(is.finite(ob)))

  # Bad eval_batch length
  sc_b_bad <- new_score(
    name = "bbad",
    order = 2L,
    eval_batch = function(subs, data, ctx = NULL, rng = NULL) rep(1, nrow(subs) + 1L)
  )
  expect_error(
    score_eval_batch(sc_b_bad, subs, data),
    "must return numeric\\(K\\) with finite entries"
  )
})

test_that("[engine] path errors when the respective R fallback is absent", {
  data <- cbind(1:5, 6:10)
  subs <- rbind(c(1L, 2L), c(3L, 4L))

  sc_only_batch <- new_score(
    name = "only_batch",
    order = 2L,
    eval_batch = function(subs, data, ctx = NULL, rng = NULL) rep(1, nrow(subs))
  )
  expect_error(
    score_eval(sc_only_batch, subs[1L, ], data),
    "No R `eval` provided; use compiled `ptr` or R `eval_batch`\\."
  )

  sc_only_eval <- new_score(
    name = "only_eval",
    order = 2L,
    eval = function(idxs, data, ctx = NULL, rng = NULL) 1
  )
  expect_error(
    score_eval_batch(sc_only_eval, subs, data),
    "No R `eval_batch` provided; use compiled `ptr` or R `eval`\\."
  )
})

test_that("[helpers] as_score wraps plain function; predicates behave", {
  f <- function(idxs, data, ctx = NULL, rng = NULL) 1.0
  sc <- as_score(f, s = 3L, name = "wrapped")
  expect_s3_class(sc, "ustat_score")
  expect_identical(score_order(sc, n = 10L, data = NULL), 3L)
  expect_true(is_score(sc))
  expect_false(is_compiled(sc))
})
