test_that("input validation catches common misuse", {
  .test_rng_scope({
    d <- .make_np_data(n = 20L, p = 2L, m = 2L)
    df <- data.frame(X1 = d$X[, 1], X2 = d$X[, 2], Y = d$y)

    expect_error(
      np_dnn(X = d$X, X0 = d$X0, x0 = d$X0[1, ], Y = d$y, k = 1L, s = 1L),
      "either `X0` or `x0`",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(X = d$X, x0 = c(NA_real_, 0), Y = d$y, k = 1L, s = 2L, se = FALSE),
      "x0.*NA",
      ignore.case = TRUE
    )
    bad_X0 <- d$X0
    bad_X0[1L, 1L] <- Inf
    expect_error(
      np_dnn(X = d$X, X0 = bad_X0, Y = d$y, k = 1L, s = 2L, se = FALSE),
      "X0.*non-finite",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(
        X = d$X, X0 = matrix(numeric(0), nrow = 0L, ncol = ncol(d$X)),
        Y = d$y, k = 1L, s = 2L, se = FALSE
      ),
      "X0.*at least 1 row",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(
        X = d$X, X0 = data.frame(X1 = "a", X2 = "b"),
        Y = d$y, k = 1L, s = 2L, se = FALSE
      ),
      "X0.*numeric",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, k = 1L, s = 1L),
      "requires `Y`",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, y = d$y, k = 1L, s = 1L),
      "`y` is only supported with `data=`",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(data = df, x_cols = c("X1", "X2"), X0 = d$X0, k = 1L, s = 1L),
      "exactly one of `y_col` or `y`",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(data = df, x_cols = c("X1", "X2"), y_col = "Y", y = d$y, X0 = d$X0, k = 1L, s = 1L),
      "exactly one of `y_col` or `y`",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = nrow(d$X), se = TRUE, jackknife = "loo"),
      "s < n",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(data = as.matrix(df), x_cols = c("X1", "X2"), y_col = "Y", X0 = d$X0),
      "data.frame",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(data = df, y_col = "Y", X0 = d$X0),
      "x_cols",
      ignore.case = TRUE
    )
    expect_error(
      np_dnn(data = df, x_cols = c("X1", "X_missing"), y_col = "Y", X0 = d$X0, k = 1L, s = 2L, se = FALSE),
      "specified by `x_cols`",
      ignore.case = TRUE
    )
    df_bad_x <- df
    df_bad_x$X1 <- as.character(df_bad_x$X1)
    expect_error(
      np_dnn(data = df_bad_x, x_cols = c("X1", "X2"), y_col = "Y", X0 = d$X0, k = 1L, s = 2L, se = FALSE),
      "X.*numeric",
      ignore.case = TRUE
    )
    # k and s validation at entry point
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 0L, s = 3L, se = FALSE),
      "positive integer"
    )
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 3L, s = 2L, se = FALSE),
      "k <= s"
    )
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1.5, s = 3L, se = FALSE),
      "integer scalar"
    )
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = "foo", s = 3L, se = FALSE),
      "integer scalar"
    )
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = 3L, se = NA),
      "se.*TRUE or FALSE"
    )
    expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = 3L, se = "yes"),
      "se.*TRUE or FALSE"
    )
    expect_error(
      np_dnn(
        X = d$X, X0 = d$X0, Y = d$y,
        method = "asymp", k = 1L, s = 3L, se = FALSE,
        normalize = NA
      ),
      "normalize.*TRUE or FALSE"
    )
    expect_error(
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

    err_na <- expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = y_na, k = 1L, s = 2L, se = FALSE),
      regexp = "NA", ignore.case = TRUE
    )
    err_inf <- expect_error(
      np_dnn(X = d$X, X0 = d$X0, Y = y_inf, k = 1L, s = 2L, se = FALSE),
      regexp = "non-finite", ignore.case = TRUE
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
    expect_error(
      np_dnn(
        X = d$X, X0 = d$X0, Y = d$y,
        distance = bad_shape, k = 1L, s = 2L, se = FALSE
      ),
      "finite numeric matrix of shape"
    )

    bad_finite <- new_distance(
      name = "np_bad_finite_distance",
      eval_batch = function(X, X0, ctx = NULL) {
        out <- matrix(0, nrow(X), nrow(X0))
        out[1L, 1L] <- Inf
        out
      }
    )
    expect_error(
      np_dnn(
        X = d$X, X0 = d$X0, Y = d$y,
        distance = bad_finite, k = 1L, s = 2L, se = FALSE
      ),
      "contains non-finite values"
    )
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
