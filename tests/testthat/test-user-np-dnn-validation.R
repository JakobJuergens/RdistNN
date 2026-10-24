.expect_np_input_error <- function(object, regexp, argument, ...) {
  err <- expect_error(object, regexp = regexp, class = "rdistnn_input_error", ...)
  expect_identical(err$code, "argument_invalid")
  expect_identical(err$argument, argument)
  expect_identical(err$operation, "np_dnn")
  invisible(err)
}

test_that("input validation catches common misuse", {
  .test_rng_scope({
    d <- .make_np_data(n = 20L, p = 2L, m = 2L)
    df <- data.frame(X1 = d$X[, 1], X2 = d$X[, 2], Y = d$y)

    .expect_np_input_error(
      np_dnn(X = d$X, Y = d$y, k = 1L, s = 2L, se = FALSE),
      "provide `X0`.*single `x0`", argument = "X0/x0"
    )
    .expect_np_input_error(
      argument = "X0/x0",
      np_dnn(X = d$X, X0 = d$X0, x0 = d$X0[1, ], Y = d$y, k = 1L, s = 1L),
      "either `X0` or `x0`",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "x0",
      np_dnn(X = d$X, x0 = c(NA_real_, 0), Y = d$y, k = 1L, s = 2L, se = FALSE),
      "x0.*NA",
      ignore.case = TRUE
    )
    bad_X0 <- d$X0
    bad_X0[1L, 1L] <- Inf
    .expect_np_input_error(
      argument = "X0",
      np_dnn(X = d$X, X0 = bad_X0, Y = d$y, k = 1L, s = 2L, se = FALSE),
      "X0.*non-finite",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "X0",
      np_dnn(
        X = d$X, X0 = matrix(numeric(0), nrow = 0L, ncol = ncol(d$X)),
        Y = d$y, k = 1L, s = 2L, se = FALSE
      ),
      "X0.*at least 1 row",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "X0",
      np_dnn(
        X = d$X, X0 = data.frame(X1 = "a", X2 = "b"),
        Y = d$y, k = 1L, s = 2L, se = FALSE
      ),
      "X0.*numeric",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "Y",
      np_dnn(X = d$X, X0 = d$X0, k = 1L, s = 1L),
      "requires `Y`",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "y",
      np_dnn(X = d$X, X0 = d$X0, y = d$y, k = 1L, s = 1L),
      "`y` is only supported with `data=`",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "y_col",
      np_dnn(data = df, x_cols = c("X1", "X2"), X0 = d$X0, k = 1L, s = 1L),
      "exactly one of `y_col` or `y`",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "y_col",
      np_dnn(data = df, x_cols = c("X1", "X2"), y_col = "Y", y = d$y, X0 = d$X0, k = 1L, s = 1L),
      "exactly one of `y_col` or `y`",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "s",
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = nrow(d$X), se = TRUE, jackknife = "loo"),
      "s < n",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "data",
      np_dnn(data = as.matrix(df), x_cols = c("X1", "X2"), y_col = "Y", X0 = d$X0),
      "data.frame",
      ignore.case = TRUE
    )
    .expect_np_input_error(
      argument = "x_cols",
      np_dnn(data = df, y_col = "Y", X0 = d$X0),
      "x_cols",
      ignore.case = TRUE
    )
    missing_column <- .expect_np_input_error(
      argument = "x_cols",
      np_dnn(data = df, x_cols = c("X1", "X_missing"), y_col = "Y", X0 = d$X0, k = 1L, s = 2L, se = FALSE),
      "specified by `x_cols`",
      ignore.case = TRUE
    )
    expect_match(missing_column$actual, "X_missing", fixed = TRUE)
    expect_match(missing_column$hint, "Check the selector.*available column names")
    df_bad_x <- df
    df_bad_x$X1 <- as.character(df_bad_x$X1)
    .expect_np_input_error(
      argument = "X",
      np_dnn(data = df_bad_x, x_cols = c("X1", "X2"), y_col = "Y", X0 = d$X0, k = 1L, s = 2L, se = FALSE),
      "X.*numeric",
      ignore.case = TRUE
    )
    # k and s validation at entry point
    .expect_np_input_error(
      argument = "k",
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 0L, s = 3L, se = FALSE),
      "positive integer"
    )
    .expect_np_input_error(
      argument = "k",
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 3L, s = 2L, se = FALSE),
      "k <= s"
    )
    .expect_np_input_error(
      argument = "k",
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1.5, s = 3L, se = FALSE),
      "integer scalar"
    )
    .expect_np_input_error(
      argument = "k",
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = "foo", s = 3L, se = FALSE),
      "integer scalar"
    )
    .expect_np_input_error(
      argument = "se",
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = 3L, se = NA),
      "se.*TRUE or FALSE"
    )
    .expect_np_input_error(
      argument = "se",
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = 3L, se = "yes"),
      "se.*TRUE or FALSE"
    )
    .expect_np_input_error(
      argument = "normalize",
      np_dnn(
        X = d$X, X0 = d$X0, Y = d$y,
        method = "asymp", k = 1L, s = 3L, se = FALSE,
        normalize = NA
      ),
      "normalize.*TRUE or FALSE"
    )
    .expect_np_input_error(
      argument = "normalize",
      np_dnn(
        X = d$X, X0 = d$X0, Y = d$y,
        method = "asymp", k = 1L, s = 3L, se = FALSE,
        normalize = "yes"
      ),
      "normalize.*TRUE or FALSE"
    )
  }, seed = 123L)
})

test_that("NA and non-finite Y values produce informative errors (matrix API)", {
  .test_rng_scope({
    d <- .make_np_data(n = 20L, p = 2L, m = 2L)
    y_na  <- d$y; y_na[5L]  <- NA_real_
    y_inf <- d$y; y_inf[5L] <- Inf

    err_na <- .expect_np_input_error(
      np_dnn(X = d$X, X0 = d$X0, Y = y_na, k = 1L, s = 2L, se = FALSE),
      regexp = "cannot contain `NA`", argument = "Y"
    )
    err_inf <- .expect_np_input_error(
      np_dnn(X = d$X, X0 = d$X0, Y = y_inf, k = 1L, s = 2L, se = FALSE),
      regexp = "cannot contain non-finite", argument = "Y"
    )
    expect_false(identical(conditionMessage(err_na), conditionMessage(err_inf)))
  }, seed = 123L)
})

test_that("np_dnn surfaces malformed custom distance outputs", {
  .test_rng_scope({
    d <- .make_np_data(n = 18L, p = 2L, m = 2L)

    bad_shape <- new_distance(
      name = "np_bad_shape_distance",
      eval_batch = function(X, X0, ctx = NULL) rep(0, nrow(X) * nrow(X0))
    )
    shape_error <- expect_error(
      np_dnn(
        X = d$X, X0 = d$X0, Y = d$y,
        distance = bad_shape, k = 1L, s = 2L, se = FALSE
      ),
      "finite numeric matrix of shape", class = "rdistnn_callback_error"
    )
    expect_identical(shape_error$code, "callback_invalid")
    expect_identical(shape_error$argument, "distance")
    expect_match(shape_error$hint, "Correct the callback.*dimensions and finite values")

    bad_finite <- new_distance(
      name = "np_bad_finite_distance",
      eval_batch = function(X, X0, ctx = NULL) {
        out <- matrix(0, nrow(X), nrow(X0))
        out[1L, 1L] <- Inf
        out
      }
    )
    finite_error <- expect_error(
      np_dnn(
        X = d$X, X0 = d$X0, Y = d$y,
        distance = bad_finite, k = 1L, s = 2L, se = FALSE
      ),
      "contains non-finite values", class = "rdistnn_callback_error"
    )
    expect_identical(finite_error$code, "callback_invalid")
    expect_identical(finite_error$argument, "distance")
    expect_match(finite_error$hint, "Correct the callback.*finite values")
  }, seed = 123L)
})

test_that("np_dnn warns on unknown return values", {
  .test_rng_scope({
    d <- .make_np_data(n = 20L, p = 2L, m = 2L)
    expect_warning(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = 2L, se = FALSE, return = c("W", "typo")),
      "Unknown `return` value"
    )
    expect_no_warning(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = 2L, se = FALSE,
             return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors"))
    )
  }, seed = 123L)
})

test_that("np_dnn rejects nonfinite signals before geometry on every input route", {
  df <- data.frame(x = 1:8, y = 2 * (1:8))
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    signal <- df$y
    signal[3L] <- bad
    bad_df <- df
    bad_df$y <- signal
    routes <- list(
      list(X = matrix(df$x), Y = signal),
      list(data = df, x_cols = "x", y = signal),
      list(data = bad_df, x_cols = "x", y_col = "y"),
      list(data = bad_df, x_cols = 1L, y_col = 2L)
    )
    for (i in seq_along(routes)) {
      probe <- .counted_l1_distance()
      err <- expect_error(do.call(np_dnn, c(routes[[i]],
        list(x0 = 0, s = 3L, se = FALSE, distance = probe$metric))),
        "cannot contain.*(NA|non-finite)", class = "rdistnn_input_error")
      expect_identical(err$code, "argument_invalid")
      expect_identical(err$argument, c("Y", "y", "y_col", "y_col")[[i]])
      expect_identical(err$operation, "np_dnn")
      expect_identical(probe$counts$setup, 0L)
      expect_identical(probe$counts$eval, 0L)
    }
  }
})

test_that("np_dnn preserves valid geometry callbacks and their failures", {
  args <- list(X = matrix(1:8), Y = 2 * (1:8), x0 = 0,
    s = 3L, se = FALSE)
  baseline <- do.call(np_dnn, c(args, list(distance = "l1")))
  probe <- .counted_l1_distance()
  actual <- do.call(np_dnn, c(args, list(distance = probe$metric)))
  expect_equal(actual$mu, baseline$mu, tolerance = 1e-14)
  expect_identical(probe$counts$setup, 1L)
  expect_identical(probe$counts$eval, 1L)

  opaque <- structure(list(message = "opaque geometry failure", call = NULL),
    class = c("opaque_geometry_error", "error", "condition"))
  for (phase in c("on_setup", "on_eval")) {
    probe <- do.call(.counted_l1_distance,
      setNames(list(function() stop(opaque)), phase))
    caught <- tryCatch(do.call(np_dnn, c(args, list(distance = probe$metric))),
      error = identity)
    expect_identical(caught, opaque)
    withr::local_options(warn = 2)
    probe <- do.call(.counted_l1_distance,
      setNames(list(function() warning("opaque geometry warning")), phase))
    promoted <- expect_error(do.call(np_dnn,
      c(args, list(distance = probe$metric))), "opaque geometry warning")
    expect_false(inherits(promoted, "rdistnn_input_error"))
  }
})
