test_that("scalar distance rankings preserve singleton training-row batch shapes", {
  setup_calls <- eval_calls <- 0L
  metric <- new_distance("scalar_distance_singleton_consumers",
    setup = function(X) {
      setup_calls <<- setup_calls + 1L
      list(scale = 2)
    },
    eval = function(X, x0, ctx) {
      eval_calls <<- eval_calls + 1L
      abs(X[, 1L] - x0[1L]) / ctx$scale
    })
  register_distance("scalar_distance_singleton_consumers", metric)
  weight <- as_rank_weight(function(ranks) rep(3, length(ranks)))
  X <- matrix(2, 1L, 1L)

  for (ranking in list(metric, "scalar_distance_singleton_consumers")) {
    for (targets in list(0, c(0, 1, 0), c(0, 0, 1, 0))) {
      setup_calls <- eval_calls <- 0L
      fit <- weights_matrix(X, X0 = matrix(targets, ncol = 1L),
        ranking = ranking, weight = weight, return = c("ranks", "ord", "dx"))
      m <- length(targets)
      expect_identical(dim(fit$W), c(1L, m))
      expect_equal(fit$W, matrix(3, 1L, m), tolerance = 0)
      expect_identical(fit$ranks, matrix(1L, 1L, m))
      expect_identical(fit$ord, matrix(1L, 1L, m))
      expect_equal(fit$dx, matrix(abs(2 - targets) / 2, 1L, m), tolerance = 0)
      expect_identical(setup_calls, 1L)
      expect_identical(eval_calls, as.integer(m))
    }
    setup_calls <- eval_calls <- 0L
    single <- weights_matrix(X, x0 = 0, ranking = ranking, weight = weight,
      return = "dx")
    expect_equal(single$W, matrix(3, 1L, 1L), tolerance = 0)
    expect_equal(single$dx, matrix(1, 1L, 1L), tolerance = 0)
    expect_identical(setup_calls, 1L)
    expect_identical(eval_calls, 1L)
  }
})

test_that("scalar distance ranking setup is shared across targets and explicit ctx is preserved", {
  setup_calls <- eval_calls <- 0L
  metric <- new_distance("scalar_distance_context_consumers",
    setup = function(X) {
      setup_calls <<- setup_calls + 1L
      list(scale = 2)
    },
    eval = function(X, x0, ctx) {
      eval_calls <<- eval_calls + 1L
      abs(X[, 1L] - x0[1L]) / ctx$scale
    })
  weight <- as_rank_weight(function(ranks, ctx) {
    multiplier <- if (is.null(ctx)) 1 else ctx$rank_scale
    multiplier * (5 - ranks)
  })
  X <- matrix(c(0, 3, 6, 9), ncol = 1L)
  X0 <- matrix(c(2, 7, 2), ncol = 1L)
  expected_ranks <- cbind(c(2L, 1L, 3L, 4L), c(4L, 3L, 1L, 2L), c(2L, 1L, 3L, 4L))
  expected_weights <- cbind(c(3, 4, 2, 1), c(1, 2, 4, 3), c(3, 4, 2, 1))
  expected_dx <- cbind(c(1, .5, 2, 3.5), c(3.5, 2, .5, 1), c(1, .5, 2, 3.5))

  automatic <- weights_matrix(X, X0, ranking = metric, weight = weight,
    return = c("ranks", "dx"))
  expect_equal(automatic$W, expected_weights, tolerance = 0)
  expect_identical(automatic$ranks, expected_ranks)
  expect_equal(automatic$dx, expected_dx, tolerance = 0)
  expect_identical(setup_calls, 1L)
  expect_identical(eval_calls, 3L)

  setup_calls <- eval_calls <- 0L
  explicit <- weights_matrix(X, X0, ranking = metric, weight = weight,
    ctx = list(scale = 4, rank_scale = 2), return = c("ranks", "dx"))
  expect_equal(explicit$W, 2 * expected_weights, tolerance = 0)
  expect_identical(explicit$ranks, expected_ranks)
  expect_equal(explicit$dx, expected_dx / 2, tolerance = 0)
  expect_identical(setup_calls, 0L)
  expect_identical(eval_calls, 3L)
})

test_that("NP and supplied CATE accept scalar distances with one training row", {
  setup_calls <- eval_calls <- 0L
  metric <- new_distance("scalar_distance_estimator_consumers",
    setup = function(X) {
      setup_calls <<- setup_calls + 1L
      list(scale = 2)
    },
    eval = function(X, x0, ctx) {
      eval_calls <<- eval_calls + 1L
      abs(X[, 1L] - x0[1L]) / ctx$scale
    })
  register_distance("scalar_distance_estimator_consumers", metric)
  X <- matrix(2, 1L, 1L)
  for (distance in list(metric, "scalar_distance_estimator_consumers")) {
    for (targets in list(0, c(0, 1, 0))) {
      X0 <- matrix(targets, ncol = 1L)
      m <- length(targets)
      for (family in c("np", "cate")) {
        setup_calls <- eval_calls <- 0L
        fit <- if (family == "np") {
          np_dnn(X = X, Y = 7, X0 = X0, distance = distance, s = 1L,
            se = FALSE, return = c("W", "dx"))
        } else {
          cate_dnn(X = X, pseudo_outcome = 7, X0 = X0, distance = distance,
            s = 1L, se = FALSE, return = c("W", "dx"))
        }
        estimate <- if (family == "np") fit$mu else fit$tau
        expect_equal(estimate, rep(7, m), tolerance = 0)
        expect_equal(fit$W, matrix(1, 1L, m), tolerance = 0)
        expect_equal(fit$dx, matrix(abs(2 - targets) / 2, 1L, m), tolerance = 0)
        expect_identical(setup_calls, 1L)
        expect_identical(eval_calls, as.integer(m))
        expect_equal(predict(fit, X0 = X0), rep(7, m), tolerance = 0)
        expect_identical(setup_calls, 1L)
        expect_identical(eval_calls, as.integer(2 * m))
      }
    }
  }
})
