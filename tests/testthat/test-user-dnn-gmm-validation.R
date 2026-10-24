# GMM input controls and malformed custom-model outputs.
test_that("gmm_dnn validates user-facing interface controls", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 70L, kx = 2L, q = 2L, seed = 28L)
  model <- .test_dr_iv_model(dgp$q)
  X <- as.matrix(dgp$data[, c("X1", "X2")])

  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = list()), "model")
  bad_model <- model
  bad_model$moment_obs <- NULL
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = bad_model), "model\\$moment_obs.*function")
  bad_model <- model
  bad_model$coef_names <- character()
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = bad_model), "model\\$coef_names")
  bad_model <- model
  bad_model$p <- 2L
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = bad_model), "model\\$p.*length\\(model\\$coef_names\\)")
  expect_error(gmm_dnn(data = dgp$data, X = X, x_cols = c("X1", "X2"),
                       x0 = dgp$x0, model = model), "data-frame interface")
  expect_error(gmm_dnn(X = X, x0 = dgp$x0, model = model), "Y")
  expect_error(gmm_dnn(X = X, model_data = dgp$data[-1L, ], x0 = dgp$x0,
                       model = model), "nrow")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = model, omega_ridge = -1), "omega_ridge")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = model, fd_step = 0), "fd_step")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = model, start = c(0, 1)), "start")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = model, start = "0"), "start.*numeric", ignore.case = TRUE)
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       X0 = matrix(dgp$x0, nrow = 1L), model = model),
               "either `X0`")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"),
                       x0 = c(NA_real_, 0), model = model, k = 1L, s = 2L),
               "x0.*NA",
               ignore.case = TRUE)
  bad_X0 <- matrix(dgp$x0, nrow = 1L)
  bad_X0[1L, 1L] <- Inf
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"),
                       X0 = bad_X0, model = model, k = 1L, s = 2L),
               "X0.*non-finite",
               ignore.case = TRUE)
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"),
                       X0 = matrix(numeric(0), nrow = 0L, ncol = 2L),
                       model = model, k = 1L, s = 2L),
               "X0.*at least 1 row",
               ignore.case = TRUE)
})

test_that("gmm_dnn validates malformed model moment outputs before target fitting", {
  withr::local_preserve_seed()
  set.seed(32L)
  n <- 35L
  df <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    Y = rnorm(n)
  )

  make_model <- function(moment_obs) {
    new_gmm_model(
      name = "moment_contract",
      setup = function(sample, n, ...) list(y = sample$Y),
      moment_obs = moment_obs,
      start = function(sample, ctx, w, x0) sum(w * ctx$y),
      coef_names = "theta",
      moment_names = "mean"
    )
  }

  fit_vector <- gmm_dnn(
    data = df,
    x_cols = c("X1", "X2"),
    x0 = c(0, 0),
    model = make_model(function(sample, theta, ctx) ctx$y - theta[1L]),
    k = 1L,
    s = 6L,
    se = FALSE
  )
  expect_equal(dim(fit_vector$theta), c(1L, 1L))
  expect_true(is.finite(fit_vector$theta[1L, 1L]))

  expect_error(
    gmm_dnn(
      data = df,
      x_cols = c("X1", "X2"),
      x0 = c(0, 0),
      model = make_model(function(sample, theta, ctx) matrix(0, nrow = nrow(sample) - 1L, ncol = 1L)),
      k = 1L,
      s = 6L,
      se = FALSE
    ),
    "model\\$moment_obs\\(\\).*one row per observation",
    ignore.case = TRUE
  )

  expect_error(
    gmm_dnn(
      data = df,
      x_cols = c("X1", "X2"),
      x0 = c(0, 0),
      model = make_model(function(sample, theta, ctx) matrix("bad", nrow = nrow(sample), ncol = 1L)),
      k = 1L,
      s = 6L,
      se = FALSE
    ),
    "model\\$moment_obs\\(\\).*numeric",
    ignore.case = TRUE
  )

  expect_error(
    gmm_dnn(
      data = df,
      x_cols = c("X1", "X2"),
      x0 = c(0, 0),
      model = make_model(function(sample, theta, ctx) {
        out <- matrix(0, nrow = nrow(sample), ncol = 1L)
        out[1L, 1L] <- NA_real_
        out
      }),
      k = 1L,
      s = 6L,
      se = FALSE
    ),
    "model\\$moment_obs\\(\\).*NA",
    ignore.case = TRUE
  )

  expect_error(
    gmm_dnn(
      data = df,
      x_cols = c("X1", "X2"),
      x0 = c(0, 0),
      model = make_model(function(sample, theta, ctx) {
        out <- matrix(0, nrow = nrow(sample), ncol = 1L)
        out[1L, 1L] <- Inf
        out
      }),
      k = 1L,
      s = 6L,
      se = FALSE
    ),
    "model\\$moment_obs\\(\\).*non-finite",
    ignore.case = TRUE
  )
})

test_that("gmm_dnn rejects malformed analytic Jacobians with callback conditions", {
  n <- 8L
  df <- data.frame(X = seq_len(n), Y = c(1, 3, 2, 5, 4, 7, 6, 8))
  model <- new_gmm_model("jacobian_contract",
    setup = function(sample, n, ...) list(y = sample$Y),
    moment_obs = function(sample, theta, ctx) ctx$y - theta[1L],
    coef_names = "theta", moment_names = "mean")
  fit_model <- function(model) gmm_dnn(data = df, x_cols = "X", x0 = 1.25,
    model = model, s = 3L, se = FALSE, jacobian = "analytic", return = "jacobian")

  failure <- tryCatch(fit_model(model), error = identity)
  expect_s3_class(failure, "rdistnn_input_error")
  expect_identical(failure$argument, "jacobian")
  expect_match(conditionMessage(failure), "analytic.*jacobian_obs")

  callbacks <- list(
    shape = function(sample, theta, ctx) matrix(-1, nrow(sample), 2L),
    type = function(sample, theta, ctx) matrix("bad", nrow(sample), 1L),
    values = function(sample, theta, ctx) c(NaN, rep(-1, nrow(sample) - 1L)))
  messages <- c(shape = "n x q matrix", type = "numeric", values = "non-finite")
  for (name in names(callbacks)) {
    model$jacobian_obs <- callbacks[[name]]
    failure <- tryCatch(fit_model(model), error = identity)
    expect_s3_class(failure, "rdistnn_callback_error")
    expect_identical(failure$code, "callback_invalid")
    expect_identical(failure$argument, "model$jacobian_obs()")
    expect_match(conditionMessage(failure), messages[[name]])
  }
})

test_that("gmm_dnn rejects invalid model-generated starts as callback failures", {
  model <- new_gmm_model("bad_start",
    setup = function(sample, n, ...) list(),
    moment_obs = function(sample, theta, ctx) rep(theta[1L], nrow(sample)),
    start = function(sample, ctx, w, x0) "0", coef_names = "theta")
  failure <- tryCatch(gmm_dnn(X = matrix(seq_len(8L)), x0 = 1.25,
    model = model, s = 3L, se = FALSE), error = identity)
  expect_s3_class(failure, "rdistnn_callback_error")
  expect_identical(failure$code, "callback_invalid")
  expect_identical(failure$argument, "model$start()")
  expect_match(conditionMessage(failure), "numeric")
})

test_that("gmm_dnn detects model moment-name mismatches", {
  withr::local_seed(123L)
  bad_model <- new_gmm_model(
    name = "bad_names",
    setup = function(sample, n, ...) list(y = sample$Y),
    moment_obs = function(sample, theta, ctx) cbind(ctx$y - theta[1L],
                                                    2 * (ctx$y - theta[1L])),
    start = function(sample, ctx, w, x0) sum(w * ctx$y),
    coef_names = "theta",
    moment_names = "only_one"
  )
  df <- data.frame(X1 = rnorm(50L), X2 = rnorm(50L), Y = rnorm(50L))

  expect_error(
    gmm_dnn(data = df, x_cols = c("X1", "X2"), x0 = c(0, 0),
            model = bad_model, k = 1L, s = 5L),
    "moment_names"
  )
})
