skip_if_not(
  all(vapply(
    c("plqr_dnn", "plqr_demo_data"),
    exists, logical(1),
    mode = "function"
  )),
  "plqr_dnn() or required helpers not found; skipping file"
)

set.seed(321)

.make_plqr_data <- function(n = 80L, p = 2L, p_d = 2L, m = 4L, tau = 0.5) {
  demo <- plqr_demo_data(n = n, p = p, p_d = p_d, tau = tau, seed = 1)
  df <- demo$data
  X <- as.matrix(df[paste0("X", seq_len(p))])
  D <- as.matrix(df[paste0("D", seq_len(p_d))])
  mD <- as.matrix(df[paste0("mD", seq_len(p_d))])
  q <- as.numeric(df$q_oracle)
  Y <- as.numeric(df$Y)
  X0 <- matrix(rnorm(m * p), m, p)
  list(
    df = df, X = X, D = D, mD = mD, q = q, Y = Y, X0 = X0,
    beta = demo$truth$beta_true, tau = tau,
    x_cols = paste0("X", seq_len(p)),
    D_cols = paste0("D", seq_len(p_d)),
    mD_cols = paste0("mD", seq_len(p_d))
  )
}

test_that("plqr_dnn matrix API matches direct g-weighted moment identity", {
  d <- .make_plqr_data()
  out <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau,
    X0 = d$X0, distance = "l2",
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = "W"
  )

  g <- (d$tau - as.numeric(d$Y <= as.vector(d$D %*% d$beta) + d$q)) * (d$D - d$mD)
  manual <- crossprod(g, out$W)

  expect_s3_class(out, "plqr_result")
  expect_equal(dim(out$moment), c(ncol(d$D), nrow(d$X0)))
  expect_equal(unname(out$moment), unname(manual), tolerance = 1e-12)
})

test_that("data-frame and matrix APIs agree for plqr_dnn", {
  d <- .make_plqr_data(n = 70L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  nms <- names(d$df)

  out_df <- plqr_dnn(
    data = d$df,
    x_cols = d$x_cols,
    D_cols = d$D_cols,
    q_col = "q_oracle",
    mD_cols = d$mD_cols,
    y_col = "Y",
    beta = d$beta,
    tau = d$tau,
    X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_df_idx <- plqr_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    D_cols = match(d$D_cols, nms),
    q_col = match("q_oracle", nms),
    mD_cols = match(d$mD_cols, nms),
    y_col = match("Y", nms),
    beta = d$beta,
    tau = d$tau,
    X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_m <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau,
    X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(out_df$moment, out_m$moment, tolerance = 1e-12)
  expect_equal(out_df$moment, out_df_idx$moment, tolerance = 1e-12)
})

test_that("plqr_dnn interfaces agree for exact method with fast jackknife", {
  d <- .make_plqr_data(n = 54L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  nms <- names(d$df)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = TRUE, jackknife = "fast"
  )
  out_df_name <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = TRUE, jackknife = "fast"
  )
  out_df_idx <- plqr_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    D_cols = match(d$D_cols, nms),
    q_col = match("q_oracle", nms),
    mD_cols = match(d$mD_cols, nms),
    y_col = match("Y", nms),
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df_name$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df_name$se, out_mat$se, tolerance = 1e-12)
  expect_equal(out_df_idx$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df_idx$se, out_mat$se, tolerance = 1e-12)
})

test_that("plqr_dnn interfaces agree for exact method with loo jackknife", {
  d <- .make_plqr_data(n = 40L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 1L, s = 2L, se = TRUE, jackknife = "loo"
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 1L, s = 2L, se = TRUE, jackknife = "loo"
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df$se, out_mat$se, tolerance = 1e-12)
})

test_that("plqr_dnn interfaces agree for asymp method with fast jackknife", {
  d <- .make_plqr_data(n = 52L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  nms <- names(d$df)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 4L, se = TRUE, jackknife = "fast"
  )
  out_df_name <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 4L, se = TRUE, jackknife = "fast"
  )
  out_df_idx <- plqr_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    D_cols = match(d$D_cols, nms),
    q_col = match("q_oracle", nms),
    mD_cols = match(d$mD_cols, nms),
    y_col = match("Y", nms),
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 4L, se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df_name$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df_name$se, out_mat$se, tolerance = 1e-12)
  expect_equal(out_df_idx$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df_idx$se, out_mat$se, tolerance = 1e-12)
})

test_that("plqr_dnn x0 path agrees across matrix and data interfaces", {
  d <- .make_plqr_data(n = 44L, p = 2L, p_d = 2L, m = 4L, tau = 0.5)
  j <- 3L

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, x0 = d$X0[j, ],
    method = "exact", k = 2L, s = 3L, se = FALSE
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, x0 = d$X0[j, ],
    method = "exact", k = 2L, s = 3L, se = FALSE
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
})

test_that("plqr_dnn return extras and meta are identical across interfaces", {
  d <- .make_plqr_data(n = 56L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 3L, s = 4L, se = FALSE,
    return = c("W", "ranks", "ord", "dx", "g", "ess", "neighbors")
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 3L, s = 4L, se = FALSE,
    return = c("W", "ranks", "ord", "dx", "g", "ess", "neighbors")
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df$W, out_mat$W, tolerance = 1e-12)
  expect_equal(out_df$ranks, out_mat$ranks)
  expect_equal(out_df$ord, out_mat$ord)
  expect_equal(out_df$dx, out_mat$dx, tolerance = 1e-12)
  expect_equal(out_df$g, out_mat$g, tolerance = 1e-12)
  expect_equal(out_df$ess, out_mat$ess, tolerance = 1e-12)
  expect_equal(out_df$neighbors, out_mat$neighbors)

  expect_identical(out_df$meta$n, out_mat$meta$n)
  expect_identical(out_df$meta$m, out_mat$meta$m)
  expect_identical(out_df$meta$p, out_mat$meta$p)
  expect_identical(out_df$meta$p_d, out_mat$meta$p_d)
  expect_identical(out_df$meta$k, out_mat$meta$k)
  expect_identical(out_df$meta$s, out_mat$meta$s)
  expect_identical(out_df$meta$method, out_mat$meta$method)
  expect_identical(out_df$meta$distance, out_mat$meta$distance)
})

test_that("plqr_dnn preserves interface equivalence with distance spec object", {
  d <- .make_plqr_data(n = 46L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  dist_spec <- distance_l2()

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    distance = dist_spec, method = "exact", k = 2L, s = 3L, se = FALSE
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, X0 = d$X0,
    distance = dist_spec, method = "exact", k = 2L, s = 3L, se = FALSE
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
})

test_that("plqr_dnn tibble input matches data.frame input", {
  skip_if_not_installed("tibble")
  d <- .make_plqr_data(n = 42L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)

  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = FALSE
  )
  out_tbl <- plqr_dnn(
    data = tibble::as_tibble(d$df), x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = FALSE
  )

  expect_equal(out_tbl$moment, out_df$moment, tolerance = 1e-12)
})

test_that("plqr_dnn interface equivalence holds in p=1, p_d=1, m=1 edge case", {
  d <- .make_plqr_data(n = 38L, p = 1L, p_d = 1L, m = 1L, tau = 0.5)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 1L, s = 2L, se = FALSE,
    return = c("W", "g")
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, q_col = "q_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 1L, s = 2L, se = FALSE,
    return = c("W", "g")
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df$W, out_mat$W, tolerance = 1e-12)
  expect_equal(out_df$g, out_mat$g, tolerance = 1e-12)
})

test_that("plqr_dnn returns finite jackknife SEs with expected shape", {
  d <- .make_plqr_data(n = 50L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)
  out <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau,
    X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_true(is.matrix(out$se))
  expect_equal(dim(out$se), dim(out$moment))
  expect_true(all(is.finite(out$se)))
  expect_true(all(out$se >= 0))
})

test_that("default jackknife mode uses fast path and matches loo reference", {
  d <- .make_plqr_data(n = 45L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)

  out_default <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = TRUE
  )
  out_fast <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_loo <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_equal(out_default$moment, out_fast$moment, tolerance = 1e-12)
  expect_equal(out_default$se, out_fast$se, tolerance = 1e-12)
  expect_equal(out_fast$se, out_loo$se, tolerance = 1e-12)
})

test_that("asymp method fast and loo jackknife paths agree", {
  d <- .make_plqr_data(n = 48L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)

  out_fast <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 3L, se = TRUE, jackknife = "fast"
  )
  out_loo <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 3L, se = TRUE, jackknife = "loo"
  )

  expect_equal(out_fast$moment, out_loo$moment, tolerance = 1e-12)
  expect_equal(out_fast$se, out_loo$se, tolerance = 1e-12)
})

test_that("return extras are populated consistently for plqr_dnn", {
  d <- .make_plqr_data(n = 55L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  out <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 3L, s = 4L, se = FALSE,
    return = c("W", "ranks", "ord", "dx", "g", "ess", "neighbors")
  )

  expect_equal(dim(out$W), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$ranks)) expect_equal(dim(out$ranks), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$ord)) expect_equal(dim(out$ord), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$dx)) expect_equal(dim(out$dx), c(nrow(d$X), nrow(d$X0)))
  expect_equal(dim(out$g), dim(d$D))
  expect_equal(unname(out$moment), unname(crossprod(out$g, out$W)), tolerance = 1e-12)
  expect_equal(out$ess, 1 / colSums(out$W^2), tolerance = 1e-12)
  expect_equal(dim(out$neighbors), c(3L, nrow(d$X0)))
  expect_true(all(out$neighbors >= 1 & out$neighbors <= nrow(d$X)))
  support_edge <- nrow(d$X) - out$meta$s + out$meta$k
  for (j in seq_len(nrow(d$X0))) {
    expect_equal(out$neighbors[, j], order(out$ranks[, j])[seq_len(out$meta$k)])
    expect_equal(sum(out$W[, j] > 0), support_edge)
    expect_gt(sum(out$W[, j] > 0), nrow(out$neighbors))
  }
})

test_that("x0 path matches corresponding single target from X0", {
  d <- .make_plqr_data(n = 44L, p = 2L, p_d = 2L, m = 4L, tau = 0.5)
  j <- 2L

  out_X0 <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_x0 <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, x0 = d$X0[j, ],
    k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(
    unname(out_x0$moment),
    unname(out_X0$moment[, j, drop = FALSE]),
    tolerance = 1e-12
  )
})

test_that("plqr_dnn validation catches shape and beta issues", {
  d <- .make_plqr_data(n = 30L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta[-1], tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "beta.*ncol\\(D\\)",
    ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD[, 1, drop = FALSE], Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "match ncol\\(D\\)|same number of columns",
    ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = nrow(d$X), se = TRUE, jackknife = "loo"
    ),
    "s < n",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = 1, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "numeric scalar in \\(0,1\\)",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0, x0 = d$X0[1, ],
      k = 1L, s = 2L, se = FALSE
    ),
    "either `X0` or `x0`",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      data = d$df,
      x_cols = d$x_cols,
      D_cols = d$D_cols,
      beta = d$beta,
      tau = d$tau,
      X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "x_cols.*D_cols.*q_col.*mD_cols.*y_col",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      data = d$df,
      x_cols = d$x_cols,
      D_cols = c(d$D_cols[1], "D_missing"),
      q_col = "q_oracle",
      mD_cols = d$mD_cols,
      y_col = "Y",
      beta = d$beta,
      tau = d$tau,
      X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "specified by `D_cols`",
    ignore.case = TRUE
  )
  df_bad_mD <- d$df
  df_bad_mD[[d$mD_cols[1]]] <- as.character(df_bad_mD[[d$mD_cols[1]]])
  expect_error(
    plqr_dnn(
      data = df_bad_mD,
      x_cols = d$x_cols,
      D_cols = d$D_cols,
      q_col = "q_oracle",
      mD_cols = d$mD_cols,
      y_col = "Y",
      beta = d$beta,
      tau = d$tau,
      X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "mD.*numeric",
    ignore.case = TRUE
  )
  # k and s validation at entry point
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 0L, s = 3L, se = FALSE
    ),
    "positive integer"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 3L, s = 2L, se = FALSE
    ),
    "k <= s"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1.5, s = 3L, se = FALSE
    ),
    "integer scalar"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = "foo", s = 3L, se = FALSE
    ),
    "integer scalar"
  )
})

test_that("plqr_dnn shape diagnostics include expected and actual values", {
  d <- .make_plqr_data(n = 28L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)
  msg <- tryCatch(
    {
      plqr_dnn(
        X = d$X, D = d$D[-1, , drop = FALSE], q = d$q, mD = d$mD, Y = d$Y,
        beta = d$beta, tau = d$tau, X0 = d$X0,
        k = 1L, s = 2L, se = FALSE
      )
      ""
    },
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, "Expected:", ignore.case = FALSE)
  expect_match(msg, "Got:", ignore.case = FALSE)
})

test_that("plqr_dnn matrix API distinguishes NA vs non-finite diagnostics", {
  d <- .make_plqr_data(n = 28L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)

  D_na <- d$D
  D_na[1, 1] <- NA_real_
  expect_error(
    plqr_dnn(
      X = d$X, D = D_na, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "D.*NA",
    ignore.case = TRUE
  )

  D_inf <- d$D
  D_inf[1, 1] <- Inf
  expect_error(
    plqr_dnn(
      X = d$X, D = D_inf, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "D.*non-finite",
    ignore.case = TRUE
  )

  q_na <- d$q
  q_na[2] <- NA_real_
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = q_na, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "q_tau\\(X\\).*NA",
    ignore.case = TRUE
  )

  mD_inf <- d$mD
  mD_inf[2, 1] <- Inf
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = mD_inf, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "mD.*non-finite",
    ignore.case = TRUE
  )

  y_inf <- d$Y
  y_inf[3] <- Inf
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = y_inf,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "Y.*non-finite",
    ignore.case = TRUE
  )
})

test_that("plqr_moment distinguishes NA vs non-finite beta diagnostics", {
  d <- .make_plqr_data(n = 30L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)

  beta_na <- d$beta
  beta_na[1] <- NA_real_
  expect_error(
    plqr_moment(d$df, beta = beta_na, tau = d$tau),
    "beta.*NA",
    ignore.case = TRUE
  )

  beta_inf <- d$beta
  beta_inf[1] <- Inf
  expect_error(
    plqr_moment(d$df, beta = beta_inf, tau = d$tau),
    "beta.*non-finite",
    ignore.case = TRUE
  )
})

test_that("plqr_dnn warns on unknown return values", {
  d <- .make_plqr_data(n = 35L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)
  expect_warning(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE, return = c("W", "typo")
    ),
    "Unknown `return` value"
  )
  expect_no_warning(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE,
      return = c("W", "ranks", "ord", "dx", "g", "ess", "neighbors")
    )
  )
})

test_that("print.plqr_result runs without error", {
  d <- .make_plqr_data(n = 35L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)
  out <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    k = 1L, s = 2L, se = FALSE
  )
  txt <- capture.output(print(out))
  expect_true(any(grepl("<plqr_result>", txt, fixed = TRUE)))
})
