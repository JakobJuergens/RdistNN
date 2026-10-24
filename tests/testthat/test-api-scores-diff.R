# Parameter-aware evaluation, analytic gradients, and finite differences.
test_that("parameter-aware evaluators and analytic gradients work", {
  y <- c(2, 3, 5, 7, 11)
  theta <- c(1.5, -0.2)
  idxs <- c(2L, 4L)
  subs <- rbind(c(1L, 2L), c(3L, 5L))

  sc <- new_score(
    name = "param_linear",
    order = 2L,
    eval = function(idxs, data, ctx, rng) sum(data[idxs]),
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
  expect_output(print(path), "<rdistnn_score_diff_path>")
  expect_equal(path$value(idxs, theta), theta[1] * sum(y[idxs]) + theta[2], tolerance = 1e-12)
  expect_equal(path$grad(idxs, theta), c(sum(y[idxs]), 1), tolerance = 1e-12)

  vb <- path$value_batch(subs, theta)
  gb <- path$grad_batch(subs, theta)
  expect_equal(vb, c(theta[1] * sum(y[c(1, 4)]) + theta[2], theta[1] * sum(y[c(2, 3)]) + theta[2]), tolerance = 1e-12)
  expect_equal(gb[, 1], c(sum(y[c(1, 4)]), sum(y[c(2, 3)])), tolerance = 1e-12)
  expect_equal(gb[, 2], c(1, 1), tolerance = 1e-12)
})

test_that("score_diff_path validates finite-difference overrides before closures are used", {
  y <- c(3, 4, 6, 9)

  sc <- new_score(
    name = "path_fd_checks",
    order = 1L,
    eval_param = function(idxs, data, theta, ctx, rng) theta[1] * data[idxs],
    grad = function(idxs, data, theta, ctx, rng) data[idxs]
  )

  expect_error(
    score_diff_path(sc, data = y, fd_step = "bad"),
    "fd_step.*positive finite numeric",
    ignore.case = TRUE
  )
  expect_error(
    score_diff_path(sc, data = y, fd_step = NA_real_),
    "fd_step.*NA",
    ignore.case = TRUE
  )
  expect_error(
    score_diff_path(sc, data = y, fd_step = c(1e-6, Inf)),
    "fd_step.*non-finite",
    ignore.case = TRUE
  )
  expect_error(
    score_diff_path(sc, data = y, fd_step = c(1e-6, 0)),
    "fd_step.*positive",
    ignore.case = TRUE
  )
  expect_error(
    score_diff_path(sc, data = y, fd_scheme = "backward"),
    "fd_scheme",
    ignore.case = TRUE
  )
  expect_error(
    score_diff_path(sc, data = y, fd_scheme = NA_character_),
    "fd_scheme",
    ignore.case = TRUE
  )

  path <- score_diff_path(
    sc,
    data = y,
    method = "fd",
    fd_step = c(1e-6, 2e-6),
    fd_scheme = "forward"
  )
  expect_equal(path$meta$fd_step, c(1e-6, 2e-6))
  expect_identical(path$meta$fd_scheme, "forward")
})

test_that("derivative helpers validate theta and finite-difference controls", {
  sc <- new_score(
    name = "theta_checks",
    order = 1L,
    eval_param = function(idxs, data, theta, ctx, rng) theta[1] * data[idxs]
  )
  y <- c(1, 2, 3)

  expect_error(score_eval_param(sc, 1L, y, numeric(0)), "non-empty numeric vector")
  expect_error(score_grad(sc, 1L, y, c(NA_real_), method = "fd"), "NA")
  expect_error(score_grad(sc, 1L, y, c(Inf), method = "fd"), "non-finite")
  expect_error(score_grad(sc, 1L, y, 1, method = "fd", fd_step = NA_real_), "cannot be `NA`")
  expect_error(score_grad(sc, 1L, y, 1, method = "fd", fd_step = Inf), "cannot be non-finite")
  expect_error(score_grad(sc, 1L, y, 1, method = "fd", fd_step = -1), "positive")
  expect_error(
    score_grad(sc, 1L, y, 1, method = "fd", fd_step = "bad"),
    "fd_step.*positive finite numeric",
    ignore.case = TRUE
  )
  expect_error(
    score_grad(sc, 1L, y, 1, method = "fd", fd_step = numeric(0)),
    "fd_step.*positive finite numeric",
    ignore.case = TRUE
  )
  expect_error(
    score_grad(sc, 1L, y, c(1, 2), method = "fd", fd_step = c(1e-6, Inf)),
    "fd_step.*non-finite",
    ignore.case = TRUE
  )
  expect_error(
    score_grad(sc, 1L, y, 1, method = "fd", fd_scheme = "backward"),
    "fd_scheme",
    ignore.case = TRUE
  )
  expect_error(
    score_grad_batch(sc, matrix(1L, nrow = 1L), y, 1, method = "fd", fd_scheme = NA_character_),
    "fd_scheme",
    ignore.case = TRUE
  )
  expect_no_error(score_grad(sc, 1L, y, c(1, 2), method = "fd", fd_step = c(1e-6)))
  expect_error(score_grad(sc, 1L, y, c(1, 2), method = "fd", fd_step = c(1e-6, -1e-6)), "positive")
})

test_that("analytic gradient outputs are validated for shape and finiteness", {
  y <- c(1, 2, 3)
  theta <- c(0.5, -1)
  subs <- rbind(c(1L, 2L), c(2L, 3L))

  sc_bad_grad_length <- new_score(
    name = "bad_grad_length",
    order = 2L,
    eval_param = function(idxs, data, theta, ctx, rng) theta[1] * sum(data[idxs]) + theta[2],
    grad = function(idxs, data, theta, ctx, rng) c(1, 2, 3)
  )
  expect_error(
    score_grad(sc_bad_grad_length, c(1L, 2L), y, theta, method = "analytic"),
    "score\\$grad.*length\\(theta\\)",
    ignore.case = TRUE
  )

  sc_bad_grad_nonfinite <- new_score(
    name = "bad_grad_nonfinite",
    order = 2L,
    eval_param = function(idxs, data, theta, ctx, rng) theta[1] * sum(data[idxs]) + theta[2],
    grad = function(idxs, data, theta, ctx, rng) c(sum(data[idxs]), Inf)
  )
  expect_error(
    score_grad(sc_bad_grad_nonfinite, c(1L, 2L), y, theta, method = "analytic"),
    "score\\$grad.*finite numeric",
    ignore.case = TRUE
  )

  sc_bad_grad_batch_shape <- new_score(
    name = "bad_grad_batch_shape",
    order = 2L,
    eval_param = function(idxs, data, theta, ctx, rng) theta[1] * sum(data[idxs]) + theta[2],
    grad_batch = function(subs, data, theta, ctx, rng) matrix(1, nrow = nrow(subs), ncol = 1L)
  )
  expect_error(
    score_grad_batch(sc_bad_grad_batch_shape, subs, y, theta, method = "analytic"),
    "score\\$grad_batch.*K x length\\(theta\\)",
    ignore.case = TRUE
  )

  sc_bad_grad_batch_nonfinite <- new_score(
    name = "bad_grad_batch_nonfinite",
    order = 2L,
    eval_param = function(idxs, data, theta, ctx, rng) theta[1] * sum(data[idxs]) + theta[2],
    grad_batch = function(subs, data, theta, ctx, rng) {
      out <- cbind(rowSums(matrix(data[subs], nrow = nrow(subs))), rep(1, nrow(subs)))
      out[2L, 1L] <- NaN
      out
    }
  )
  expect_error(
    score_grad_batch(sc_bad_grad_batch_nonfinite, subs, y, theta, method = "analytic"),
    "score\\$grad_batch.*finite matrix",
    ignore.case = TRUE
  )
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
