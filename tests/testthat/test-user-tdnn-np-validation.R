test_that("input validation catches common misuse", {
  .test_rng_scope({
    d  <- .make_tdnn_data(n = 25L, p = 2L, m = 2L)
    df <- data.frame(X1 = d$X[, 1], X2 = d$X[, 2], Y = d$y)

    # Mutually exclusive X0/x0
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, x0 = d$X0[1L, ], Y = d$y, s1 = 1L, s2 = 2L),
      "either `X0` or `x0`",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(X = d$X, x0 = c(NA_real_, 0), Y = d$y, s1 = 1L, s2 = 2L, se = FALSE),
      "x0.*NA",
      ignore.case = TRUE
    )
    bad_X0 <- d$X0
    bad_X0[1L, 1L] <- Inf
    expect_error(
      np_tdnn(X = d$X, X0 = bad_X0, Y = d$y, s1 = 1L, s2 = 2L, se = FALSE),
      "X0.*non-finite",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(
        X = d$X, X0 = matrix(numeric(0), nrow = 0L, ncol = ncol(d$X)),
        Y = d$y, s1 = 1L, s2 = 2L, se = FALSE
      ),
      "X0.*at least 1 row",
      ignore.case = TRUE
    )
    # Matrix API requires Y
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, s1 = 1L, s2 = 2L),
      "requires `Y`|`Y` is required",
      ignore.case = TRUE
    )
    # y only supported in data= path
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, y = d$y, s1 = 1L, s2 = 2L),
      "`y` is only supported with `data=`",
      ignore.case = TRUE
    )
    # data= requires exactly one of y_col/y
    expect_error(
      np_tdnn(data = df, x_cols = c("X1", "X2"), X0 = d$X0, s1 = 1L, s2 = 2L),
      "exactly one of `y_col` or `y`",
      ignore.case = TRUE
    )
    # s1 >= s2
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 3L, s2 = 2L, se = FALSE),
      "s1 < s2",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 3L, s2 = 3L, se = FALSE),
      "s1 < s2",
      ignore.case = TRUE
    )
    # s1 must be a positive integer
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 0L, s2 = 2L, se = FALSE),
      "positive integer",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 1.5, s2 = 3L, se = FALSE),
      "integer scalar",
      ignore.case = TRUE
    )
    # Logical controls must be explicit scalar TRUE/FALSE values.
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 1L, s2 = 2L, se = NA),
      "se.*TRUE or FALSE",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 1L, s2 = 2L, se = "yes"),
      "se.*TRUE or FALSE",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 1L, s2 = 2L, se = c(TRUE, FALSE)),
      "se.*TRUE or FALSE",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(
        X = d$X, X0 = d$X0, Y = d$y,
        s1 = 1L, s2 = 2L, method = "asymp", se = FALSE, normalize = NA
      ),
      "normalize.*TRUE or FALSE",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(
        X = d$X, X0 = d$X0, Y = d$y,
        s1 = 1L, s2 = 2L, method = "asymp", se = FALSE, normalize = "yes"
      ),
      "normalize.*TRUE or FALSE",
      ignore.case = TRUE
    )
    expect_error(
      np_tdnn(
        X = d$X, X0 = d$X0, Y = d$y,
        s1 = 1L, s2 = 2L, method = "asymp", se = FALSE, normalize = c(TRUE, FALSE)
      ),
      "normalize.*TRUE or FALSE",
      ignore.case = TRUE
    )
    # Jackknife feasibility: s2 >= n
    expect_error(
      np_tdnn(
        X = d$X, X0 = d$X0, Y = d$y,
        s1 = 1L, s2 = nrow(d$X), se = TRUE, jackknife = "loo"
      ),
      "s2 < n",
      ignore.case = TRUE
    )
    # data= requires x_cols
    expect_error(
      np_tdnn(data = df, y_col = "Y", X0 = d$X0, s1 = 1L, s2 = 2L),
      "x_cols",
      ignore.case = TRUE
    )
    # data= must be a data.frame
    expect_error(
      np_tdnn(
        data = as.matrix(df), x_cols = c("X1", "X2"), y_col = "Y",
        X0 = d$X0, s1 = 1L, s2 = 2L
      ),
      "data.frame",
      ignore.case = TRUE
    )
  }, seed = 123L)
})

test_that("np_tdnn warns on unknown return values", {
  .test_rng_scope({
    d <- .make_tdnn_data(n = 25L, m = 2L)
    expect_warning(
      np_tdnn(
        X = d$X, X0 = d$X0, Y = d$y,
        s1 = 1L, s2 = 2L, se = FALSE, return = c("W", "typo")
      ),
      "Unknown `return` value"
    )
    expect_no_warning(
      np_tdnn(
        X = d$X, X0 = d$X0, Y = d$y,
        s1 = 1L, s2 = 2L, se = FALSE,
        return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")
      )
    )
  }, seed = 123L)
})

test_that("np_tdnn jackknife feasibility diagnostics include expected and actual values", {
  .test_rng_scope({
    d <- .make_tdnn_data(n = 25L, p = 2L, m = 2L)
    msg <- tryCatch(
      {
        np_tdnn(
          X = d$X, X0 = d$X0, Y = d$y,
          s1 = 1L, s2 = nrow(d$X), se = TRUE, jackknife = "loo"
        )
        ""
      },
      error = function(e) conditionMessage(e)
    )
    expect_match(msg, "Expected:", ignore.case = FALSE)
    expect_match(msg, "Got:", ignore.case = FALSE)
  }, seed = 123L)
})

test_that("np_tdnn rejects nonfinite signals before geometry on every input route", {
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
      err <- expect_error(do.call(np_tdnn, c(routes[[i]],
        list(x0 = 0, s1 = 2L, s2 = 3L, se = FALSE, distance = probe$metric))),
        "cannot contain.*(NA|non-finite)", class = "rdistnn_input_error")
      expect_identical(err$code, "argument_invalid")
      expect_identical(err$argument, c("Y", "y", "y_col", "y_col")[[i]])
      expect_identical(err$operation, "np_tdnn")
      expect_identical(probe$counts$setup, 0L)
      expect_identical(probe$counts$eval, 0L)
    }
  }
})

test_that("np_tdnn preserves valid geometry callbacks and their failures", {
  args <- list(X = matrix(1:8), Y = 2 * (1:8), x0 = 0,
    s1 = 2L, s2 = 3L, se = FALSE)
  baseline <- do.call(np_tdnn, c(args, list(distance = "l1")))
  probe <- .counted_l1_distance()
  actual <- do.call(np_tdnn, c(args, list(distance = probe$metric)))
  expect_equal(actual$mu, baseline$mu, tolerance = 1e-14)
  expect_identical(probe$counts$setup, 1L)
  expect_identical(probe$counts$eval, 1L)

  opaque <- structure(list(message = "opaque geometry failure", call = NULL),
    class = c("opaque_geometry_error", "error", "condition"))
  for (phase in c("on_setup", "on_eval")) {
    probe <- do.call(.counted_l1_distance,
      setNames(list(function() stop(opaque)), phase))
    caught <- tryCatch(do.call(np_tdnn, c(args, list(distance = probe$metric))),
      error = identity)
    expect_identical(caught, opaque)
    withr::local_options(warn = 2)
    probe <- do.call(.counted_l1_distance,
      setNames(list(function() warning("opaque geometry warning")), phase))
    promoted <- expect_error(do.call(np_tdnn,
      c(args, list(distance = probe$metric))), "opaque geometry warning")
    expect_false(inherits(promoted, "rdistnn_input_error"))
  }
})
