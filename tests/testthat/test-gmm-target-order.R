.audit_gmm_square_model <- function(start = NULL, shifted = FALSE) {
  new_gmm_model(
    name = "target_order_square",
    setup = function(sample, n, ...) {
      list(center = if (shifted) sample$center else rep(0, n))
    },
    moment_obs = function(sample, theta, ctx) {
      matrix((theta[1L] - ctx$center)^2 - 4, ncol = 1L)
    },
    jacobian_obs = function(sample, theta, ctx) {
      matrix(2 * (theta[1L] - ctx$center), ncol = 1L)
    },
    start = start,
    coef_names = "theta",
    moment_names = "square"
  )
}

test_that("GMM calls model starts for each target independently of target order", {
  calls <- numeric()
  model <- .audit_gmm_square_model(start = function(sample, ctx, w, x0) {
    calls <<- c(calls, x0[1L])
    if (x0[1L] < 0) -1 else 1
  })
  fit <- function(targets) {
    gmm_dnn(X = matrix(seq(-2, 2, length.out = 8L)),
      X0 = matrix(targets), model = model, s = 3L,
      se = FALSE, return = character())
  }
  targets <- c(-1, 1, -1)
  batched <- fit(targets)
  expect_identical(calls, targets)
  expect_equal(as.numeric(coef(batched)), c(-2, 2, -2), tolerance = 1e-5)
  expect_true(all(batched$convergence$convergence_first == 0L))
  expect_true(all(batched$convergence$convergence_second == 0L))

  scalar <- vapply(targets, function(target) as.numeric(coef(fit(target))), numeric(1L))
  expect_equal(as.numeric(coef(batched)), scalar, tolerance = 1e-10)
  perm <- c(2L, 3L, 1L)
  expect_equal(as.numeric(coef(fit(targets[perm]))), scalar[perm], tolerance = 1e-10)
})

test_that("GMM caller starts override model callbacks at every target", {
  model <- .audit_gmm_square_model(start = function(...) {
    stop("The model start must not run when caller starts are supplied.")
  })
  fit <- function(targets, start) {
    gmm_dnn(X = matrix(seq(-2, 2, length.out = 8L)),
      X0 = matrix(targets), model = model, s = 3L, start = start,
      se = FALSE, return = character())
  }
  targets <- c(-1, 1, -1)
  expect_equal(as.numeric(coef(fit(targets, 1))), rep(2, 3L), tolerance = 1e-5)
  starts <- matrix(c(-1, 1, -1), nrow = 1L)
  expected <- c(-2, 2, -2)
  expect_equal(as.numeric(coef(fit(targets, starts))), expected, tolerance = 1e-5)
  perm <- c(2L, 3L, 1L)
  expect_equal(as.numeric(coef(fit(targets[perm], starts[, perm, drop = FALSE]))),
    expected[perm], tolerance = 1e-5)
})

test_that("GMM default starts preserve scalar and permuted target results", {
  model <- .audit_gmm_square_model(shifted = TRUE)
  fit <- function(targets) {
    gmm_dnn(X = matrix(1:8), X0 = matrix(targets),
      model_data = data.frame(center = c(3, rep(.5, 7L))),
      model = model, s = 8L, se = FALSE, return = character())
  }
  targets <- c(1, 8, 1)
  scalar <- vapply(targets, function(target) as.numeric(coef(fit(target))), numeric(1L))
  expect_equal(scalar, c(1, -1.5, 1), tolerance = 1e-5)
  expect_equal(as.numeric(coef(fit(targets))), scalar, tolerance = 1e-10)
  perm <- c(2L, 3L, 1L)
  expect_equal(as.numeric(coef(fit(targets[perm]))), scalar[perm], tolerance = 1e-10)
})

.audit_gmm_positive_domain_model <- function(start = NULL) {
  new_gmm_model("positive_domain",
    setup = function(sample, n, ...) list(y = sample$y),
    moment_obs = function(sample, theta, ctx) {
      matrix(log(theta[1L]) - ctx$y, ncol = 1L)
    },
    jacobian_obs = function(sample, theta, ctx) {
      matrix(1 / theta[1L], length(ctx$y), 1L)
    },
    start = start, coef_names = "positive_parameter")
}

test_that("GMM probes positive-domain moments at explicit target starts", {
  .with_distance_tie_warning({
  model <- .audit_gmm_positive_domain_model(start = function(...) {
    stop("Explicit starts must bypass the model callback.")
  })
  targets <- c(2, 6, 2)
  y <- seq(.1, .8, length.out = 8L)
  for (start in list(1, matrix(c(1, 2, 1), nrow = 1L))) {
    fit <- gmm_dnn(X = matrix(1:8), X0 = matrix(targets),
      model_data = data.frame(y = y), model = model,
      start = start, s = 3L, se = FALSE)
    expect_equal(as.numeric(coef(fit)), as.numeric(exp(crossprod(y, fit$W))), tolerance = 1e-5)
    expect_true(all(fit$convergence$convergence_second == 0L))
  }
  })
})

test_that("GMM reuses the first callback start after its dimension probe", {
  .with_distance_tie_warning({
  calls <- numeric()
  model <- .audit_gmm_positive_domain_model(start = function(sample, ctx, w, x0) {
    calls <<- c(calls, x0[1L])
    1
  })
  targets <- c(6, 2, 6)
  y <- seq(.1, .8, length.out = 8L)
  fit <- gmm_dnn(X = matrix(1:8), X0 = matrix(targets),
    model_data = data.frame(y = y), model = model, s = 3L, se = FALSE)
  expect_identical(calls, targets)
  expect_equal(as.numeric(coef(fit)), as.numeric(exp(crossprod(y, fit$W))), tolerance = 1e-5)
  })
})

test_that("GMM keeps zero fallback and propagates target start failures", {
  .with_distance_tie_warning({
  evaluations <- numeric()
  model <- new_gmm_model("start_failure_control",
    setup = function(sample, n, ...) list(y = sample$y),
    moment_obs = function(sample, theta, ctx) {
      evaluations <<- c(evaluations, theta[1L])
      matrix(theta[1L] - ctx$y, ncol = 1L)
    }, coef_names = "mean")
  fit <- gmm_dnn(X = matrix(1:8), X0 = matrix(c(2, 6)),
    model_data = data.frame(y = 1:8), model = model, s = 3L, se = FALSE)
  expect_identical(evaluations[1L], 0)
  expect_true(all(is.finite(coef(fit))))

  calls <- numeric()
  model$start <- function(sample, ctx, w, x0) {
    calls <<- c(calls, x0[1L])
    if (x0[1L] == 2) stop("first target start failed")
    1
  }
  expect_error(gmm_dnn(X = matrix(1:8), X0 = matrix(c(2, 6)),
    model_data = data.frame(y = 1:8), model = model, s = 3L, se = FALSE),
    "first target start failed")
  expect_identical(calls, 2)
  })
})

test_that("GMM propagates later callback failures after a valid positive-domain probe", {
  .with_distance_tie_warning({
  calls <- numeric()
  model <- .audit_gmm_positive_domain_model(start = function(sample, ctx, w, x0) {
    calls <<- c(calls, x0[1L])
    if (x0[1L] == 2) stop("first target start failed")
    1
  })
  targets <- c(6, 2, 7)
  y <- seq(.1, .8, length.out = 8L)
  expect_error(gmm_dnn(X = matrix(1:8), X0 = matrix(targets),
    model_data = data.frame(y = y), model = model, s = 3L, se = FALSE),
    "first target start failed")
  expect_identical(calls, c(6, 2))
  })
})

test_that("GMM rejects invalid-domain explicit starts without trying later targets", {
  .with_distance_tie_warning({
  model <- .audit_gmm_positive_domain_model(start = function(...) {
    stop("Explicit starts must bypass the model callback.")
  })
  y <- seq(.1, .8, length.out = 8L)
  failure <- tryCatch(gmm_dnn(X = matrix(1:8), X0 = matrix(c(2, 6)),
    model_data = data.frame(y = y), model = model,
    start = matrix(c(0, 1), nrow = 1L), s = 3L, se = FALSE), error = identity)
  expect_s3_class(failure, "rdistnn_callback_error")
  expect_identical(failure$code, "callback_invalid")
  expect_identical(failure$argument, "model$moment_obs()")
  expect_match(conditionMessage(failure), "non-finite")
  })
})

test_that("GMM does not probe zero after the first callback start fails", {
  calls <- evaluations <- numeric()
  model <- new_gmm_model("all_start_failures",
    setup = function(sample, n, ...) list(y = sample$y),
    moment_obs = function(sample, theta, ctx) {
      evaluations <<- c(evaluations, theta[1L])
      matrix(theta[1L] - ctx$y, ncol = 1L)
    },
    start = function(sample, ctx, w, x0) {
      calls <<- c(calls, x0[1L])
      stop("target start failed")
    }, coef_names = "mean")
  targets <- c(2.25, 6.25, 7.25)
  expect_error(gmm_dnn(X = matrix(1:8), X0 = matrix(targets),
    model_data = data.frame(y = 1:8), model = model, s = 3L, se = FALSE),
    "target start failed")
  expect_identical(calls, targets[1L])
  expect_length(evaluations, 0L)
})
