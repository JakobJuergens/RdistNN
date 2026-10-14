# tests/testthat/test-api-scores.R
# ------------------------------------------------------------------------------

test_that("constructors and predicates: new_score / as_score / is_score / is_compiled / print", {
  # minimal eval-only score
  sc1 <- new_score(
    name = "sc1",
    order = 2L,
    eval = function(idxs, data, ctx, rng) sum(data[idxs])
  )
  expect_true(is_score(sc1))
  expect_false(is_compiled(sc1))
  expect_output(print(sc1), "<ustat_score: sc1>")
  expect_output(print(sc1), "order: 2")
  expect_output(print(sc1), "R eval ")

  # as_score wrapper
  sc2 <- as_score(function(idxs, data, ctx, rng) prod(data[idxs]), s = 3L, name = "sc2")
  expect_true(is_score(sc2))
  expect_false(is_compiled(sc2))
  expect_output(print(sc2), "<ustat_score: sc2>")

  # as_score wrapper also accepts short signatures (idxs, data)
  sc2_short <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "sc2_short")
  expect_equal(score_eval(sc2_short, c(1L, 3L), 1:5, NULL), 4)
  subs <- matrix(c(1L, 2L, 2L, 5L), nrow = 2L, byrow = TRUE)
  expect_equal(score_eval_batch(sc2_short, subs, 1:5, NULL), c(3, 7))

  # as_score wrapper accepts variadic signatures and forwards only supported args
  sc2_dots <- as_score(function(idxs, data, ...) sum(data[idxs]) + 1, s = 2L, name = "sc2_dots")
  expect_equal(score_eval(sc2_dots, c(2L, 4L), 1:5, NULL), 7)

  # pretend "compiled"
  sc3 <- new_score(
    name = "sc3",
    order = 2L,
    ptr = TRUE # any non-NULL marks compiled
  )
  expect_true(is_compiled(sc3))
  expect_output(print(sc3), "compiled\\(ptr\\)")

  # warning for non-symmetric
  expect_warning(
    new_score(name = "sc4", order = 2L, eval = function(i, d, c, r) 0, symmetric = FALSE),
    "Score marked non-symmetric; engine assumes symmetry\\."
  )
})

test_that("setup and dynamic order resolution work", {
  calls <- list(setup = 0L)

  sc <- new_score(
    name = "dyn",
    order = function(n, data, ctx) if (n >= 5) 3L else 2L,
    setup = function(data, ...) {
      calls$setup <<- calls$setup + 1L
      list(ok = TRUE)
    },
    eval = function(idxs, data, ctx, rng) sum(data[idxs])
  )

  data <- 1:6
  ctx <- score_setup(sc, data)
  expect_equal(calls$setup, 1L)
  expect_true(is.list(ctx))
  expect_equal(score_order(sc, length(data), data, ctx), 3L)
  expect_equal(score_order(sc, 4L, data, ctx), 2L)
})

test_that("score constructors and dynamic order reject malformed slots early", {
  expect_error(
    new_score(name = "bad_order_zero", order = 0L, eval = function(idxs, data, ctx, rng) 0),
    "order.*positive",
    ignore.case = TRUE
  )
  expect_error(
    new_score(name = c("a", "b"), order = 1L, eval = function(idxs, data, ctx, rng) 0),
    "name.*non-empty character scalar",
    ignore.case = TRUE
  )
  expect_error(
    new_score(name = "bad_setup", order = 1L, setup = 1, eval = function(idxs, data, ctx, rng) 0),
    "setup.*function",
    ignore.case = TRUE
  )
  expect_error(
    new_score(name = "bad_eval", order = 1L, eval = 1),
    "eval.*function",
    ignore.case = TRUE
  )
  expect_error(
    new_score(name = "bad_eval_batch", order = 1L, eval_batch = 1),
    "eval_batch.*function",
    ignore.case = TRUE
  )

  y <- 1:5
  sc_dyn_fractional <- new_score(
    name = "dyn_fractional",
    order = function(n, data, ctx) 2.5,
    eval = function(idxs, data, ctx, rng) 0
  )
  expect_error(
    score_order(sc_dyn_fractional, length(y), y, NULL),
    "score\\$order\\(\\).*integer scalar"
  )

  sc_dyn_na <- new_score(
    name = "dyn_na",
    order = function(n, data, ctx) NA_integer_,
    eval = function(idxs, data, ctx, rng) 0
  )
  expect_error(
    score_eval(sc_dyn_na, c(1L, 2L), y, NULL),
    "score\\$order\\(\\).*integer scalar"
  )
  expect_error(score_order(sc_dyn_na, 0L, y, NULL), "n.*positive", ignore.case = TRUE)
})

test_that("eval prefers eval_batch when available; batch-only score routes single subset through batch", {
  set.seed(1)
  n <- 8L
  s <- 3L
  y <- rnorm(n)

  calls <- list(batch = 0L, scalar = 0L)

  sc_batch_only <- new_score(
    name = "batch_only",
    order = s,
    eval_batch = function(subs, data, ctx, rng) {
      calls$batch <<- calls$batch + 1L
      # sum of y over each subset (row)
      rowSums(matrix(data[subs], nrow = nrow(subs), ncol = ncol(subs)))
    }
  )

  # Single subset → should use eval_batch path
  idxs <- c(1L, 4L, 5L)
  val1 <- score_eval(sc_batch_only, idxs, y, NULL)
  expect_equal(val1, sum(y[idxs]))
  expect_equal(calls$batch, 1L)

  # Many subsets
  subs <- rbind(c(1L, 4L, 5L), c(2L, 3L, 8L), c(6L, 7L, 5L))
  val2 <- score_eval_batch(sc_batch_only, subs, y, NULL)
  expect_equal(val2, c(sum(y[c(1, 4, 5)]), sum(y[c(2, 3, 8)]), sum(y[c(6, 7, 5)])))
  expect_equal(calls$batch, 2L)
})

test_that("batch fallback uses scalar eval row-wise when eval_batch is missing", {
  set.seed(2)
  n <- 7L
  s <- 2L
  y <- rnorm(n)

  calls <- list(batch = 0L, scalar = 0L)

  sc_eval_only <- new_score(
    name = "eval_only",
    order = s,
    eval = function(idxs, data, ctx, rng) {
      calls$scalar <<- calls$scalar + 1L
      sum(data[idxs]^2)
    }
  )

  subs <- rbind(c(1L, 2L), c(3L, 4L), c(6L, 7L))
  out <- score_eval_batch(sc_eval_only, subs, y, NULL)
  expect_equal(out, c(sum(y[1:2]^2), sum(y[3:4]^2), sum(y[6:7]^2)))
  expect_equal(calls$scalar, nrow(subs))

  # Single subset via score_eval still works (uses eval directly here)
  expect_equal(score_eval(sc_eval_only, c(5L, 7L), y, NULL), sum(y[c(5, 7)]^2))
})

test_that("subset validation in batch: shape, bounds, and empties are caught", {
  n <- 6L
  s <- 3L
  y <- rep(1, n)
  sc <- new_score(name = "v", order = s, eval = function(i, d, c, r) sum(d[i]))

  # wrong number of columns
  expect_error(
    score_eval_batch(sc, matrix(1L, nrow = 2L, ncol = s - 1L), y, NULL),
    "must have 3 columns"
  )

  # out-of-bounds indices
  expect_error(
    score_eval_batch(sc, matrix(c(1L, 2L, n + 1L), 1L, s, byrow = TRUE), y, NULL),
    "indices must lie in 1..n"
  )

  # empty
  expect_error(
    score_eval_batch(sc, matrix(integer(), 0L, s), y, NULL),
    "must have at least one row"
  )

  # NA vs non-finite are distinct diagnostics
  expect_error(
    score_eval_batch(sc, matrix(c(1, 2, NA_real_), 1L, s, byrow = TRUE), y, NULL),
    "subs.*NA",
    ignore.case = TRUE
  )
  expect_error(
    score_eval_batch(sc, matrix(c(1, 2, Inf), 1L, s, byrow = TRUE), y, NULL),
    "subs.*non-finite",
    ignore.case = TRUE
  )
})

test_that("score evaluators validate subsets before calling custom evaluators", {
  y <- 1:5

  scalar_calls <- 0L
  sc_scalar_ignores_indices <- new_score(
    name = "scalar_ignores_indices",
    order = 2L,
    eval = function(idxs, data, ctx, rng) {
      scalar_calls <<- scalar_calls + 1L
      0
    }
  )

  expect_error(
    score_eval(sc_scalar_ignores_indices, c(1L, 6L), y, NULL),
    "idxs.*1\\.\\.n",
    ignore.case = TRUE
  )
  expect_error(
    score_eval(sc_scalar_ignores_indices, c(1, NA_real_), y, NULL),
    "idxs.*NA",
    ignore.case = TRUE
  )
  expect_error(
    score_eval(sc_scalar_ignores_indices, c(1, 2.5), y, NULL),
    "idxs.*integer indices",
    ignore.case = TRUE
  )
  expect_identical(scalar_calls, 0L)

  expect_no_error(score_eval(sc_scalar_ignores_indices, c(2L, 2L), y, NULL))
  expect_identical(scalar_calls, 1L)

  batch_calls <- 0L
  sc_batch_ignores_indices <- new_score(
    name = "batch_ignores_indices",
    order = 2L,
    eval_batch = function(subs, data, ctx, rng) {
      batch_calls <<- batch_calls + 1L
      rep(0, nrow(subs))
    }
  )

  expect_error(
    score_eval_batch(sc_batch_ignores_indices, matrix(c(1L, 6L), nrow = 1L), y, NULL),
    "subs.*1\\.\\.n",
    ignore.case = TRUE
  )
  expect_error(
    score_eval_batch(sc_batch_ignores_indices, matrix(c(1, Inf), nrow = 1L), y, NULL),
    "subs.*non-finite",
    ignore.case = TRUE
  )
  expect_error(
    score_eval_batch(sc_batch_ignores_indices, matrix(c(1, 2.5), nrow = 1L), y, NULL),
    "subs.*integer indices",
    ignore.case = TRUE
  )
  expect_identical(batch_calls, 0L)

  valid_subs <- rbind(c(1L, 1L), c(2L, 5L))
  expect_equal(score_eval_batch(sc_batch_ignores_indices, valid_subs, y, NULL), c(0, 0))
  expect_identical(batch_calls, 1L)
})

test_that("registry: register_score(name, score), get_score(), and list_scores()", {
  scA <- new_score(name = "scA", order = 2L, eval = function(i, d, c, r) 1)

  nmA <- register_score("scA", scA)
  expect_identical(nmA, "scA")
  gotA <- get_score("scA")
  expect_true(is_score(gotA))
  expect_identical(gotA$name, "scA")

  # (name, score): register under a different key; internal name stays as constructed
  scB <- new_score(name = "scB_inner", order = 2L, eval = function(i, d, c, r) 2)
  nmB <- register_score("altB", scB)
  expect_identical(nmB, "altB")
  gotB <- get_score("altB")
  expect_true(is_score(gotB))
  # The registry key is "altB", but the score's own name remains "scB_inner"
  expect_identical(gotB$name, "scB_inner")

  # listing
  nn <- list_scores()
  expect_true("scA" %in% nn)
  expect_true("altB" %in% nn)
})

test_that("get_score diagnostics include argument and expected/got details", {
  expect_error(get_score(1), "name", ignore.case = TRUE)
  msg <- tryCatch(
    {
      get_score("__missing_score__")
      NA_character_
    },
    error = function(e) e$message
  )
  expect_match(msg, "No score registered", ignore.case = TRUE)
  expect_match(msg, "Expected:", fixed = TRUE)
  expect_match(msg, "Got:", fixed = TRUE)
})

test_that("errors when neither eval_batch nor eval are present (and no ptr)", {
  # Construct with a dummy ptr to pass validation, then strip ptr to simulate a broken object.
  sc_bad <- new_score(name = "bad", order = 2L, ptr = TRUE)
  sc_bad$ptr <- NULL

  expect_error(
    score_eval(sc_bad, c(1L, 2L), 1:5, NULL),
    "No R `eval` or `eval_batch` provided"
  )
  expect_error(
    score_eval_batch(sc_bad, matrix(c(1L, 2L), 1L, 2L), 1:5, NULL),
    "No R `eval_batch` or `eval` provided"
  )
})

test_that("return shape and finiteness checks are enforced", {
  y <- 1:5
  sc_bad_len <- new_score(
    name = "bad_len",
    order = 2L,
    eval_batch = function(subs, data, ctx, rng) rep(0, nrow(subs) + 1L)
  )
  expect_error(
    score_eval_batch(sc_bad_len, matrix(c(1L, 2L), 1L, 2L), y, NULL),
    "must return numeric\\(K\\)"
  )

  sc_bad_inf <- new_score(
    name = "bad_inf",
    order = 2L,
    eval = function(idxs, data, ctx, rng) Inf
  )
  expect_error(
    score_eval(sc_bad_inf, c(1L, 2L), y, NULL),
    "must return a finite scalar"
  )
})

test_that("randomized flag is carried through (print only smoke test)", {
  sc <- new_score(name = "rand", order = 2L, eval = function(i, d, c, r) 0, randomized = TRUE)
  expect_output(print(sc), "\\(randomized\\)")
})
