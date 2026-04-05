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
})

test_that("registry: register_score() single-arg and (name, score), get_score(), list_scores(), score_names()", {
  scA <- new_score(name = "scA", order = 2L, eval = function(i, d, c, r) 1)

  # single-arg; uses scA$name
  nmA <- register_score(scA)
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
  expect_true(all(score_names() %in% list_scores()))
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

test_that("parameter-aware evaluators and analytic gradients work", {
  y <- c(2, 3, 5, 7, 11)
  theta <- c(1.5, -0.2)
  idxs <- c(2L, 4L)
  subs <- rbind(c(1L, 2L), c(3L, 5L))

  sc <- new_score(
    name = "param_linear",
    order = 2L,
    eval = function(idxs, data, ctx, rng) sum(data[idxs]), # legacy path still present
    eval_param = function(idxs, data, theta, ctx, rng) {
      theta[1] * sum(data[idxs]) + theta[2]
    },
    eval_param_batch = function(subs, data, theta, ctx, rng) {
      theta[1] * rowSums(matrix(data[subs], nrow = nrow(subs), ncol = ncol(subs))) + theta[2]
    },
    grad = function(idxs, data, theta, ctx, rng) c(sum(data[idxs]), 1),
    grad_batch = function(subs, data, theta, ctx, rng) {
      s1 <- rowSums(matrix(data[subs], nrow = nrow(subs), ncol = ncol(subs)))
      cbind(s1, rep(1, length(s1)))
    }
  )

  expect_equal(score_eval_param(sc, idxs, y, theta), theta[1] * sum(y[idxs]) + theta[2], tolerance = 1e-12)
  expect_equal(
    score_eval_param_batch(sc, subs, y, theta),
    theta[1] * c(sum(y[c(1, 2)]), sum(y[c(3, 5)])) + theta[2],
    tolerance = 1e-12
  )

  expect_equal(score_grad(sc, idxs, y, theta, method = "analytic"), c(sum(y[idxs]), 1), tolerance = 1e-12)
  G <- score_grad_batch(sc, subs, y, theta, method = "analytic")
  expect_equal(G[, 1], c(sum(y[c(1, 2)]), sum(y[c(3, 5)])), tolerance = 1e-12)
  expect_equal(G[, 2], c(1, 1), tolerance = 1e-12)
})

test_that("finite-difference gradient path works when analytic gradients are absent", {
  y <- c(1.2, -0.4, 0.8, 2.0, -1.5)
  theta <- c(0.7, -0.3)
  idxs <- c(1L, 3L, 4L)
  subs <- rbind(c(1L, 3L, 4L), c(2L, 4L, 5L))

  sc_fd <- new_score(
    name = "param_nonlinear_fd",
    order = 3L,
    eval_param = function(idxs, data, theta, ctx, rng) {
      s1 <- sum(data[idxs])
      m1 <- mean(data[idxs])
      (theta[1]^2) * s1 + sin(theta[2]) * m1
    },
    fd_step = 1e-6,
    fd_scheme = "central"
  )

  g_ref <- c(
    2 * theta[1] * sum(y[idxs]),
    cos(theta[2]) * mean(y[idxs])
  )
  g_hat <- score_grad(sc_fd, idxs, y, theta, method = "fd")
  expect_equal(g_hat, g_ref, tolerance = 1e-5)

  G_hat <- score_grad_batch(sc_fd, subs, y, theta, method = "fd")
  G_ref <- rbind(
    c(2 * theta[1] * sum(y[subs[1, ]]), cos(theta[2]) * mean(y[subs[1, ]])),
    c(2 * theta[1] * sum(y[subs[2, ]]), cos(theta[2]) * mean(y[subs[2, ]]))
  )
  expect_equal(G_hat, G_ref, tolerance = 1e-5)

  expect_error(score_grad(sc_fd, idxs, y, theta, method = "analytic"), "No analytic derivative")
})

test_that("score_diff_path exposes modular value/gradient closures", {
  y <- c(3, 4, 6, 9)
  theta <- c(2.0, 0.5)
  idxs <- c(1L, 4L)
  subs <- rbind(c(1L, 4L), c(2L, 3L))

  sc <- new_score(
    name = "path_score",
    order = 2L,
    eval_param = function(idxs, data, theta, ctx, rng) theta[1] * sum(data[idxs]) + theta[2],
    grad = function(idxs, data, theta, ctx, rng) c(sum(data[idxs]), 1)
  )

  path <- score_diff_path(sc, data = y, method = "auto")
  expect_output(print(path), "<score_diff_path>")
  expect_equal(path$value(idxs, theta), theta[1] * sum(y[idxs]) + theta[2], tolerance = 1e-12)
  expect_equal(path$grad(idxs, theta), c(sum(y[idxs]), 1), tolerance = 1e-12)

  vb <- path$value_batch(subs, theta)
  gb <- path$grad_batch(subs, theta)
  expect_equal(vb, c(theta[1] * sum(y[c(1, 4)]) + theta[2], theta[1] * sum(y[c(2, 3)]) + theta[2]), tolerance = 1e-12)
  expect_equal(gb[, 1], c(sum(y[c(1, 4)]), sum(y[c(2, 3)])), tolerance = 1e-12)
  expect_equal(gb[, 2], c(1, 1), tolerance = 1e-12)
})

test_that("derivative helpers validate theta and finite-difference controls", {
  sc <- new_score(
    name = "theta_checks",
    order = 1L,
    eval_param = function(idxs, data, theta, ctx, rng) theta[1] * data[idxs]
  )
  y <- c(1, 2, 3)

  expect_error(score_eval_param(sc, 1L, y, numeric(0)), "non-empty numeric vector")
  expect_error(score_grad(sc, 1L, y, c(NA_real_), method = "fd"), "finite")
  expect_error(score_grad(sc, 1L, y, 1, method = "fd", fd_step = -1), "positive")
  expect_no_error(score_grad(sc, 1L, y, c(1, 2), method = "fd", fd_step = c(1e-6)))
  expect_error(score_grad(sc, 1L, y, c(1, 2), method = "fd", fd_step = c(1e-6, -1e-6)), "positive")
})

test_that("finite differences can use eval/eval_batch when they accept theta", {
  y <- c(2, 5, 8, 10)
  theta <- c(1.3, -0.4)
  idxs <- c(1L, 3L)
  subs <- rbind(c(1L, 3L), c(2L, 4L))

  sc <- new_score(
    name = "theta_in_eval",
    order = 2L,
    eval = function(idxs, data, theta, ctx, rng) {
      theta[1] * sum(data[idxs]) + theta[2]
    },
    eval_batch = function(subs, data, theta, ctx, rng) {
      theta[1] * rowSums(matrix(data[subs], nrow = nrow(subs), ncol = ncol(subs))) + theta[2]
    }
  )

  expect_equal(score_eval_param(sc, idxs, y, theta), theta[1] * sum(y[idxs]) + theta[2], tolerance = 1e-12)
  expect_equal(
    score_eval_param_batch(sc, subs, y, theta),
    theta[1] * c(sum(y[c(1, 3)]), sum(y[c(2, 4)])) + theta[2],
    tolerance = 1e-12
  )

  g <- score_grad(sc, idxs, y, theta, method = "fd")
  expect_equal(g, c(sum(y[idxs]), 1), tolerance = 1e-5)
})

test_that("parameter-aware path errors if eval/eval_batch do not declare theta", {
  y <- c(1, 2, 3, 4)
  idxs <- c(1L, 2L)
  theta <- c(0.5)

  sc <- new_score(
    name = "theta_missing",
    order = 2L,
    eval = function(idxs, data, ctx, rng) sum(data[idxs]),
    eval_batch = function(subs, data, ctx, rng) rowSums(matrix(data[subs], nrow = nrow(subs)))
  )

  expect_error(
    score_eval_param(sc, idxs, y, theta),
    regexp = "does not declare `theta`/`param`"
  )
  expect_error(
    score_eval_param_batch(sc, rbind(idxs), y, theta),
    regexp = "does not declare `theta`/`param`"
  )
})
