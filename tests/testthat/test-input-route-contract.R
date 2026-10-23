test_that("core input routes reject ignored direct inputs and selectors", {
  d <- data.frame(x = 1:8, y = 1:8)
  for (fun in list(np_dnn, np_tdnn)) {
    order_args <- if (identical(fun, np_dnn)) list(s = 3L) else list(s1 = 2L, s2 = 4L)
    args <- c(list(data = d, x_cols = "x", y_col = "y", x0 = 0, se = FALSE), order_args)
    expect_error(do.call(fun, c(args, list(X = matrix(1:8)))), "X")
    expect_error(do.call(fun, c(args, list(Y = 1:8))), "Y")
    expect_error(do.call(fun, c(list(X = matrix(1:8), Y = 1:8, x_cols = "x", x0 = 0,
      se = FALSE), order_args)), "data")
    args$y_col <- NULL; args$y <- d$y
    expect_s3_class(do.call(fun, args), "np_result")
  }
  expect_error(cate_dnn(data = d, x_cols = "x", X = matrix(1:8),
    pseudo_outcome = 1:8, x0 = 0, s = 3, se = FALSE), "X")
  expect_s3_class(cate_dnn(data = d, x_cols = "x", pseudo_outcome = 1:8,
    x0 = 0, s = 3, se = FALSE), "cate_result")
})

test_that("GMM selectors require data while separate matrix model data stays supported", {
  model <- new_gmm_model("input_route_mean",
    setup = function(sample, n, ...) list(g = sample$g),
    moment_obs = function(sample, theta, ctx) matrix(ctx$g - theta[1L], ncol = 1L),
    jacobian_obs = function(sample, theta, ctx) matrix(-1, length(ctx$g), 1L),
    start = function(sample, ctx, w, x0) sum(w * ctx$g),
    coef_names = "theta")
  d <- data.frame(x = 1:8, g = c(4, 1, 3, 7, 2, 5, 8, 6))
  args <- list(X = matrix(d$x), model_data = d["g"], model = model,
    x0 = 0, s = 3L, se = FALSE)
  expect_error(do.call(gmm_dnn, c(args, list(x_cols = "not_a_column"))),
    "data.*column selectors")
  matrix_fit <- do.call(gmm_dnn, args)
  frame_fit <- gmm_dnn(data = d, x_cols = "x", model = model, x0 = 0, s = 3L, se = FALSE)
  expect_equal(matrix_fit$theta, frame_fit$theta, tolerance = 1e-12)
})

test_that("tabular selector range errors precede integer conversion", {
  data <- data.frame(x = 1:4, y = c(2, 1, 4, 3))
  for (selector in c(2^31, -2^31, 1e100, -1e100,
                     .Machine$double.xmax, -.Machine$double.xmax)) {
    for (argument in c("x_cols", "y_col")) {
      args <- list(data = data, x_cols = 1L, y_col = 2L,
        x0 = 0, s = 2L, se = FALSE)
      args[[argument]] <- selector
      warnings <- character()
      error <- withCallingHandlers(tryCatch({
        do.call(np_dnn, args)
        NULL
      }, error = identity), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      expect_s3_class(error, "error")
      expect_match(if (inherits(error, "error")) conditionMessage(error) else "",
        paste0(argument, ".*out-of-range"))
      expect_length(warnings, 0L)
    }
  }
  for (selector in c(1L, ncol(data))) {
    tabular <- np_dnn(data = data, x_cols = selector, y_col = 2L,
      x0 = 0, s = 2L, se = FALSE)
    explicit <- np_dnn(X = as.matrix(data[selector]), Y = data$y,
      x0 = 0, s = 2L, se = FALSE)
    expect_equal(tabular$mu, explicit$mu, tolerance = 0)
    expect_identical(tabular$.model$x_cols, names(data)[selector])
  }
})
