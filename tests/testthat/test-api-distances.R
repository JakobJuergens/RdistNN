# tests/testthat/test-distance-api.R
# ------------------------------------------------------------------------------

test_that("as_distance() wraps a function and is_distance() recognizes it", {
  withr::local_preserve_seed()
  f <- function(X, x0, ctx) {
    rowSums((X - matrix(x0, nrow(X), ncol(X), byrow = TRUE))^2)
  }
  spec <- as_distance(f, name = "quad")
  expect_true(is_distance(spec))
  expect_s3_class(spec, "rdistnn_distance")

  # single eval
  set.seed(1)
  X <- matrix(rnorm(12), 4, 3)
  x0 <- c(0, 0, 0)
  d <- distance_eval(spec, X, x0)
  expect_type(d, "double")
  expect_length(d, nrow(X))

  # batch via fallback (spec has only eval)
  X0 <- matrix(c(0, 0, 0, 1, 1, 1), 2, 3, byrow = TRUE)
  D <- distance_eval_batch(spec, X, X0)
  expect_equal(dim(D), c(nrow(X), nrow(X0)))
  expect_equal(D[, 1], d)
})

test_that("new_distance() builds validated specs", {
  spec <- new_distance(
    name = "new_quad",
    eval = function(X, x0, ctx = NULL) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
  )
  expect_true(is_distance(spec))
  expect_equal(spec$name, "new_quad")

  expect_error(
    new_distance(name = "bad_new"),
    "must provide at least one of `eval` or `eval_batch`"
  )
})

test_that("new_distance() rejects malformed names and optional slots", {
  ok_eval <- function(X, x0, ctx = NULL) rep(0, nrow(X))
  ok_batch <- function(X, X0, ctx = NULL) matrix(0, nrow(X), nrow(X0))

  expect_error(
    new_distance(name = c("a", "b"), eval = ok_eval),
    "name.*non-empty character scalar",
    ignore.case = TRUE
  )
  expect_error(
    new_distance(name = NA_character_, eval = ok_eval),
    "name.*non-empty character scalar",
    ignore.case = TRUE
  )
  expect_error(
    new_distance(name = 1, eval = ok_eval),
    "name.*non-empty character scalar",
    ignore.case = TRUE
  )
  expect_error(
    new_distance(name = "bad_setup", eval = ok_eval, setup = 1),
    "setup.*function",
    ignore.case = TRUE
  )
  expect_error(
    new_distance(name = "bad_eval", eval = 1, eval_batch = ok_batch),
    "eval.*function",
    ignore.case = TRUE
  )
  expect_error(
    new_distance(name = "bad_eval_batch", eval = ok_eval, eval_batch = 1),
    "eval_batch.*function",
    ignore.case = TRUE
  )
})

test_that("register_distance()/get_distance()/distance_names() work", {
  withr::local_seed(123L)
  # Register a simple absolute-distance spec under a fresh name
  nm <- paste0("abs_", sample.int(1e6, 1))
  spec <- as_distance(function(X, x0, ctx) rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE))), name = nm)

  invisible(register_distance(nm, spec))
  got <- get_distance(nm)
  expect_true(is_distance(got))
  expect_equal(got$name, nm)
  expect_true(nm %in% distance_names())

  # Also allow registering a function directly
  nm2 <- paste0("lin_", sample.int(1e6, 1))
  fun <- function(X, x0, ctx) drop(X %*% abs(x0)) # silly linear form
  invisible(register_distance(nm2, fun))
  got2 <- get_distance(nm2)
  expect_true(is_distance(got2))
  expect_equal(got2$name, nm2)

  # Unknown name errors
  expect_error(get_distance("__definitely_missing__"), class = "rdistnn_input_error")
  msg <- tryCatch(
    {
      get_distance("__definitely_missing__")
      NA_character_
    },
    error = function(e) e$message
  )
  expect_match(msg, "Expected:", fixed = TRUE)
  expect_match(msg, "Got:", fixed = TRUE)
})

test_that("as_distance() accepts user functions with or without ctx argument", {
  X <- matrix(1:12, nrow = 4)
  x0 <- c(0, 0, 0)

  spec_ctx <- as_distance(function(X, x0, ctx) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2))
  spec_no_ctx <- as_distance(function(X, x0) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2))

  d1 <- distance_eval(spec_ctx, X, x0)
  d2 <- distance_eval(spec_no_ctx, X, x0)
  expect_equal(d1, d2)
})

test_that("distance_eval_batch() calls setup() exactly once (eval-only fallback)", {
  withr::local_preserve_seed()
  # spec with setup counter; only eval provided (so batch uses fallback loop)
  env <- new.env(parent = emptyenv())
  env$count <- 0L

  spec <- list(
    name = "with_setup_eval_only",
    setup = function(X) {
      env$count <- env$count + 1L
      list(scale = 2)
    },
    eval = function(X, x0, ctx) {
      # use ctx$scale to prove setup ran
      rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2) * ctx$scale
    }
  )
  class(spec) <- "rdistnn_distance"

  set.seed(2)
  X <- matrix(rnorm(20), 10, 2)
  X0 <- matrix(rnorm(6), 3, 2)

  D <- distance_eval_batch(spec, X, X0, ctx = NULL)
  expect_equal(dim(D), c(nrow(X), nrow(X0)))
  expect_equal(env$count, 1L) # setup called once at batch level
})

test_that("distance_eval() falls back to eval_batch when only eval_batch exists", {
  spec <- list(
    name = "batch_only",
    eval_batch = function(X, X0, ctx) {
      # squared L2, returns n x m
      vapply(seq_len(nrow(X0)), function(j) {
        rowSums((X - matrix(X0[j, ], nrow = nrow(X), ncol = ncol(X), byrow = TRUE))^2)
      }, numeric(nrow(X)))
    }
  )
  class(spec) <- "rdistnn_distance"

  X <- matrix(1:12, 4, 3)
  x0 <- c(0, 1, 2)

  d <- distance_eval(spec, X, x0) # falls back to eval_batch internally
  expect_length(d, nrow(X))

  D <- distance_eval_batch(spec, X, matrix(x0, nrow = 1L), ctx = NULL)
  expect_equal(as.numeric(D[, 1L]), d)
})

test_that("distance_eval honors provided ctx (no auto-setup)", {
  cnt <- new.env(parent = emptyenv())
  cnt$setup <- 0L
  spec <- list(
    name = "ctx_override",
    setup = function(X) {
      cnt$setup <- cnt$setup + 1L
      list(scale = 99)
    },
    eval = function(X, x0, ctx) {
      stopifnot(is.list(ctx), !is.null(ctx$scale))
      rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2) * ctx$scale
    }
  )
  class(spec) <- "rdistnn_distance"

  X <- matrix(1, 5, 2)
  x0 <- c(0, 0)

  # user-provided ctx prevents setup()
  ctx <- list(scale = 3)
  d <- distance_eval(spec, X, x0, ctx = ctx)
  expect_true(all(d == 3 * rowSums((X - 0)^2)))
  expect_equal(cnt$setup, 0L)
})

test_that("distance_eval/Batch basic shape checks and errors", {
  # wrong-length outputs from user fun -> rejected by as_distance()
  bad <- as_distance(function(X, x0, ctx) 1, name = "bad_len")
  X <- matrix(0, 3, 2)
  x0 <- c(0, 0)
  expect_error(distance_eval(bad, X, x0), "must return a finite numeric vector of length")

  # mismatched ncol in batch triggers error
  ok <- as_distance(function(X, x0, ctx) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2),
    name = "ok"
  )
  X0 <- matrix(0, 2, 3) # wrong p
  expect_error(distance_eval_batch(ok, X, X0), "ncol\\(X0\\) must match ncol\\(X\\)")
})

test_that("distance evaluators validate targets before calling custom functions", {
  X <- matrix(0, 3, 2)

  eval_calls <- 0L
  spec_eval_ignores_target <- new_distance(
    name = "ignore_target_eval",
    eval = function(X, x0, ctx = NULL) {
      eval_calls <<- eval_calls + 1L
      rep(0, nrow(X))
    }
  )

  expect_error(
    distance_eval(spec_eval_ignores_target, X, c(NA_real_, 0)),
    "x0.*NA",
    ignore.case = TRUE
  )
  expect_error(
    distance_eval(spec_eval_ignores_target, X, c(Inf, 0)),
    "x0.*non-finite",
    ignore.case = TRUE
  )
  expect_error(
    distance_eval(spec_eval_ignores_target, X, c("a", "b")),
    "x0.*numeric",
    ignore.case = TRUE
  )
  expect_identical(eval_calls, 0L)

  expect_equal(distance_eval(spec_eval_ignores_target, X, c(0, 0)), rep(0, nrow(X)))
  expect_identical(eval_calls, 1L)

  batch_calls <- 0L
  spec_batch_ignores_target <- new_distance(
    name = "ignore_target_batch",
    eval_batch = function(X, X0, ctx = NULL) {
      batch_calls <<- batch_calls + 1L
      matrix(0, nrow(X), nrow(X0))
    }
  )

  expect_error(
    distance_eval_batch(spec_batch_ignores_target, X, matrix(c(NA_real_, 0), nrow = 1L)),
    "X0.*NA",
    ignore.case = TRUE
  )
  expect_error(
    distance_eval_batch(spec_batch_ignores_target, X, matrix(c(Inf, 0), nrow = 1L)),
    "X0.*non-finite",
    ignore.case = TRUE
  )
  expect_error(
    distance_eval_batch(spec_batch_ignores_target, X, matrix(c("a", "b"), nrow = 1L)),
    "X0.*numeric",
    ignore.case = TRUE
  )
  expect_identical(batch_calls, 0L)

  D <- distance_eval_batch(spec_batch_ignores_target, X, matrix(c(0, 0, 1, 1), nrow = 2L, byrow = TRUE))
  expect_equal(D, matrix(0, nrow(X), 2L))
  expect_identical(batch_calls, 1L)
})

test_that("registry specs actually evaluate after registration", {
  withr::local_seed(123L)
  nm <- paste0("toy_", sample.int(1e6, 1))
  fun <- function(X, x0, ctx) rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  register_distance(nm, fun)

  spec <- get_distance(nm)
  X <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3), 4, 2, byrow = TRUE)
  x0 <- c(1, 1)
  d <- distance_eval(spec, X, x0)
  expect_equal(d, c(2, 0, 2, 4))
})

test_that("distance_eval errors on x0 length mismatch", {
  spec <- as_distance(function(X, x0, ctx) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2),
    name = "quad"
  )
  X <- matrix(0, 4, 3)
  x0 <- c(1, 2) # wrong p
  expect_error(distance_eval(spec, X, x0), "x0 length must equal ncol\\(X\\)")
})

test_that("distance_eval rejects non-finite outputs", {
  spec <- as_distance(function(X, x0, ctx) rep(NaN, nrow(X)), name = "bad_nonfinite")
  X <- matrix(0, 3, 2)
  x0 <- c(0, 0)
  expect_error(distance_eval(spec, X, x0), "contains non-finite values")
})

test_that("distance_eval_batch respects provided ctx (no setup) in eval-only fallback", {
  env <- new.env(parent = emptyenv())
  env$count <- 0L
  spec <- list(
    name = "eval_only_ctx_passthrough",
    setup = function(X) {
      env$count <- env$count + 1L
      list(scale = 9)
    },
    eval = function(X, x0, ctx) {
      stopifnot(is.list(ctx), !is.null(ctx$scale))
      rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2) * ctx$scale
    }
  )
  class(spec) <- "rdistnn_distance"

  X <- matrix(1:12, 4, 3)
  X0 <- matrix(c(0, 0, 0, 1, 1, 1), 2, 3, byrow = TRUE)

  dctx <- list(scale = 2)
  D <- distance_eval_batch(spec, X, X0, ctx = dctx)
  expect_equal(dim(D), c(4L, 2L))
  expect_true(all(D[, 1] == rowSums((X - 0)^2) * 2))
  expect_equal(env$count, 0L) # setup not called because ctx was provided
})

test_that("distance_names() returns sorted registry names", {
  withr::local_seed(123L)
  nmA <- paste0("zA_", sample.int(1e6, 1))
  nmB <- paste0("aB_", sample.int(1e6, 1))
  register_distance(nmA, function(X, x0, ctx) rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE))))
  register_distance(nmB, function(X, x0, ctx) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2))
  dn <- distance_names()
  expect_true(all(sort(dn) == dn))
  expect_true(all(c(nmA, nmB) %in% dn))
})

test_that("distance eval_batch path validates shape and finiteness", {
  X <- matrix(1:12, nrow = 4)
  X0 <- matrix(c(0, 0, 0, 1, 1, 1), nrow = 2, byrow = TRUE)

  bad_shape <- list(
    name = "bad_shape",
    eval_batch = function(X, X0, ctx) rep(0, nrow(X) * nrow(X0))
  )
  class(bad_shape) <- "rdistnn_distance"
  expect_error(
    distance_eval_batch(bad_shape, X, X0),
    "must return a finite numeric matrix of shape"
  )

  bad_finite <- list(
    name = "bad_finite",
    eval_batch = function(X, X0, ctx) {
      out <- matrix(0, nrow(X), nrow(X0))
      out[1, 1] <- Inf
      out
    }
  )
  class(bad_finite) <- "rdistnn_distance"
  expect_error(
    distance_eval_batch(bad_finite, X, X0),
    "must return a finite numeric matrix of shape"
  )
})

test_that("register_distance validates malformed specs and named list specs", {
  withr::local_seed(123L)
  # malformed: no eval and no eval_batch
  expect_error(
    register_distance("bad_no_eval", list(name = "bad_no_eval")),
    "must provide at least one of `eval` or `eval_batch`"
  )

  # malformed name in two-arg form
  expect_error(
    register_distance("", list(eval = function(X, x0, ctx) rowSums(X))),
    "non-empty character scalar"
  )

  expect_error(
    register_distance(
      "bad_inner_name",
      list(name = c("a", "b"), eval = function(X, x0, ctx = NULL) rowSums(X))
    ),
    "name.*non-empty character scalar",
    ignore.case = TRUE
  )

  # named list spec in two-arg form is accepted and normalized
  nm <- paste0("list_spec_", sample.int(1e6, 1))
  spec <- list(
    eval = function(X, x0, ctx = NULL) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
  )
  register_distance(nm, spec)
  got <- get_distance(nm)
  expect_true(inherits(got, "rdistnn_distance"))
  expect_equal(got$name, nm)
})

test_that("distance_eval validates raw list specs before evaluation", {
  X <- matrix(0, 3, 2)
  raw_bad_setup <- list(
    name = "raw_bad_setup",
    setup = 1,
    eval = function(X, x0, ctx = NULL) rep(0, nrow(X))
  )

  expect_true(is_distance(raw_bad_setup))
  expect_error(
    distance_eval(raw_bad_setup, X, c(0, 0)),
    "setup.*function",
    ignore.case = TRUE
  )
})

test_that("distance_eval_batch errors when eval looks like batch but eval_batch is absent", {
  spec <- list(
    name = "ambiguous_eval",
    eval = function(X, X0, ctx = NULL) matrix(0, nrow(X), nrow(X0))
  )
  class(spec) <- "rdistnn_distance"

  X <- matrix(1:12, nrow = 4)
  X0 <- matrix(c(0, 0, 0, 1, 1, 1), nrow = 2, byrow = TRUE)
  expect_error(distance_eval_batch(spec, X, X0), "No R `eval_batch` or `eval` provided")
})

test_that("distance callback return contracts and thrown failures stay distinct", {
  X <- matrix(1:4, ncol = 1)
  bad <- new_distance("bad_result", eval = function(X, x0, ctx) "bad")
  error <- tryCatch(distance_eval(bad, X, 0), error = identity)
  expect_s3_class(error, "rdistnn_callback_error")
  expect_identical(error$code, "callback_invalid")
  expect_identical(error$argument, "distance")
  sentinel <- structure(list(message = "callback sentinel", call = NULL),
    class = c("callback_sentinel", "error", "condition"))
  thrown <- new_distance("thrown", eval = function(X, x0, ctx) stop(sentinel))
  expect_identical(tryCatch(distance_eval(thrown, X, 0), error = identity), sentinel)
})

test_that("distance lookup errors preserve the calling argument", {
  X <- matrix(1:4, ncol = 1)
  calls <- list(name = function() get_distance("l2_scaled_typo"),
    distance = function() np_dnn(X = X, Y = 1:4, x0 = 0, s = 2,
      distance = "l2_scaled_typo", se = FALSE))
  for (argument in names(calls)) {
    error <- tryCatch(calls[[argument]](), error = identity)
    expect_s3_class(error, "rdistnn_input_error")
    expect_identical(error$argument, argument)
    expect_match(conditionMessage(error), "registered distance", fixed = TRUE)
  }
})


test_that("weight wrappers retain their public distance argument in lookup errors", {
  for (fun in list(weights_dnn_exact, weights_dnn_asymp)) {
    error <- tryCatch(fun(matrix(1:4), x0 = 0, k = 1, s = 2,
      distance = "l2_scaled_typo"), error = identity)
    expect_s3_class(error, "rdistnn_input_error")
    expect_identical(error$argument, "distance")
  }
})
