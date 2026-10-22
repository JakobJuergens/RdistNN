.retained_geometry_gmm_mean_model <- function() {
  new_gmm_model(
    name = "retained_geometry_mean",
    setup = function(sample, n, ...) list(y = sample$Y),
    moment_obs = function(sample, theta, ctx) matrix(ctx$y - theta[1L], ncol = 1L),
    jacobian_obs = function(sample, theta, ctx) matrix(-1, length(ctx$y), 1L),
    start = function(sample, ctx, w, x0) sum(w * ctx$y),
    coef_names = "theta"
  )
}

test_that("GMM weight reporting retains its resolved metric and fitted context", {
  original <- get_distance("l2")
  on.exit(register_distance("l2", original), add = TRUE)
  setup_calls <- 0L
  metric <- new_distance("gmm_retained_geometry",
    setup = function(X) {
      setup_calls <<- setup_calls + 1L
      list(scale = max(X))
    },
    eval = function(X, x0, ctx) abs(X[, 1L] - x0[1L]) / ctx$scale
  )
  register_distance("l2", metric)
  fit <- gmm_dnn(X = matrix(1:8), X0 = matrix(c(3.5, 6.5)),
    model_data = data.frame(Y = c(1, 3, 2, 7, 4, 6, 5, 8)),
    model = .retained_geometry_gmm_mean_model(), distance = "l2", s = 3L,
    return = "W")
  expect_identical(setup_calls, 1L)
  before <- plot_data_weights(fit)
  expect_identical(setup_calls, 1L)
  expect_identical(before$obs_id[before$target_id == 1L], c(3L, 4L, 2L, 5L, 1L, 6L, 7L, 8L))
  register_distance("l2", new_distance("replacement",
    eval = function(X, x0, ctx = NULL) 10 - abs(X[, 1L] - x0[1L])))
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(plot_data_weights(restored), before)
  expect_identical(plot_data_weights(fit), before)
  expect_identical(setup_calls, 1L)
})

test_that("GMM explicit metric contexts remain frozen during weight reporting", {
  metric <- new_distance("gmm_explicit_geometry",
    setup = function(X) stop("An explicit metric context must bypass setup."),
    eval = function(X, x0, ctx) {
      rowSums(abs(sweep(X, 2L, x0, "-")) *
        matrix(ctx$inv_scale, nrow(X), ncol(X), byrow = TRUE))
    }
  )
  X <- cbind(1:8, c(3, 1, 7, 2, 4, 8, 6, 5))
  target <- matrix(c(3, 4), nrow = 1L)
  fit <- gmm_dnn(X = X, X0 = target, ctx = list(inv_scale = c(1, .25)),
    model_data = data.frame(Y = c(1, 3, 2, 7, 4, 6, 5, 8)),
    model = .retained_geometry_gmm_mean_model(), distance = metric, s = 3L,
    return = "W")
  expected <- order(abs(X[, 1L] - 3) + abs(X[, 2L] - 4) / 4)
  actual <- plot_data_weights(unserialize(serialize(fit, NULL)))
  expect_identical(actual$obs_id, expected)
  expect_identical(actual$weight, as.numeric(fit$W[expected, 1L]))
})
