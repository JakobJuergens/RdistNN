# tests/testthat/test-user-cate-dnn-data-interface.R

skip_if_not(
  all(vapply(
    c("cate_dnn", "cate_demo_data"),
    exists, logical(1),
    mode = "function"
  )),
  "cate_dnn() or demo-data helpers not found; skipping file"
)

test_that("data-frame AIPW path: psi matches manual formula", {
  d <- .make_df_data()
  nms <- names(d$df)
  out <- cate_dnn(
    data = d$df,
    x_cols = d$x_cols,
    y_col = "Y",
    t_col = "D",
    m0_col = "m0_oracle",
    m1_col = "m1_oracle",
    e_col = "e_oracle",
    X0 = d$X0,
    k = 3L, s = 5L, method = "exact", se = FALSE,
    return = c("W", "psi")
  )
  out_idx <- cate_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    y_col = match("Y", nms),
    t_col = match("D", nms),
    m0_col = match("m0_oracle", nms),
    m1_col = match("m1_oracle", nms),
    e_col = match("e_oracle", nms),
    X0 = d$X0,
    k = 3L, s = 5L, method = "exact", se = FALSE,
    return = c("W", "psi")
  )

  expect_equal(out$psi, d$psi_aipw, tolerance = 1e-12)
  expect_equal(out_idx$psi, out$psi, tolerance = 1e-12)
  tau_manual <- as.numeric(crossprod(out$psi, out$W))
  expect_equal(out$tau, tau_manual, tolerance = 1e-12)
  expect_length(out$tau, nrow(d$X0))
})

test_that("data-frame Robinson path: psi matches manual formula", {
  d <- .make_df_data()
  out <- cate_dnn(
    data = d$df,
    x_cols = d$x_cols,
    y_col = "Y",
    t_col = "D",
    m_col = "m_oracle",
    e_col = "e_oracle",
    X0 = d$X0,
    k = 3L, s = 5L, method = "exact", se = FALSE,
    return = c("W", "psi")
  )

  expect_equal(out$psi, d$psi_robinson, tolerance = 1e-12)
  expect_equal(out$psi, as.numeric(d$df$rho_oracle), tolerance = 1e-12)
  tau_manual <- as.numeric(crossprod(out$psi, out$W))
  expect_equal(out$tau, tau_manual, tolerance = 1e-12)
})

test_that("AIPW and Robinson paths produce different psi (and tau) on the same data", {
  d <- .make_df_data(seed = 11L)
  run <- function(use_aipw) {
    if (use_aipw) {
      cate_dnn(
        data = d$df, x_cols = d$x_cols,
        y_col = "Y", t_col = "D",
        m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
        X0 = d$X0, k = 2L, s = 4L, se = FALSE, return = "psi"
      )
    } else {
      cate_dnn(
        data = d$df, x_cols = d$x_cols,
        y_col = "Y", t_col = "D",
        m_col = "m_oracle", e_col = "e_oracle",
        X0 = d$X0, k = 2L, s = 4L, se = FALSE, return = "psi"
      )
    }
  }
  out_a <- run(TRUE)
  out_r <- run(FALSE)
  expect_true(all(is.finite(out_a$psi)))
  expect_true(all(is.finite(out_r$psi)))
  expect_false(isTRUE(all.equal(out_a$psi, out_r$psi)))
})

test_that("nuisance return slot: AIPW path returns m0/m1/e; Robinson path returns m/e", {
  d <- .make_df_data(seed = 13L)

  out_a <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
    X0 = d$X0, k = 2L, s = 3L, se = FALSE,
    return = "nuisance"
  )
  expect_true(!is.null(out_a$nuisance))
  expect_setequal(names(out_a$nuisance), c("m0", "m1", "e"))
  expect_equal(out_a$nuisance$m0, as.numeric(d$df$m0_oracle), tolerance = 1e-12)
  expect_equal(out_a$nuisance$m1, as.numeric(d$df$m1_oracle), tolerance = 1e-12)
  expect_equal(out_a$nuisance$e, as.numeric(d$df$e_oracle), tolerance = 1e-12)

  out_r <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m_col = "m_oracle", e_col = "e_oracle",
    X0 = d$X0, k = 2L, s = 3L, se = FALSE,
    return = "nuisance"
  )
  expect_true(!is.null(out_r$nuisance))
  expect_setequal(names(out_r$nuisance), c("m", "e"))
  expect_equal(out_r$nuisance$m, as.numeric(d$df$m_oracle), tolerance = 1e-12)
  expect_equal(out_r$nuisance$e, as.numeric(d$df$e_oracle), tolerance = 1e-12)
})

test_that("ESS equals 1/colSums(W^2) exactly", {
  set.seed(99)
  n <- 30
  p <- 2
  m <- 4
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  psi <- rnorm(n)
  out <- cate_dnn(
    X = X, X0 = X0, pseudo_outcome = psi,
    distance = "l2", k = 2L, s = 4L, method = "exact", se = FALSE,
    return = c("W", "ess")
  )
  expect_equal(out$ess, 1 / colSums(out$W^2), tolerance = 1e-12)
})

test_that("data-frame API errors: conflicting / missing outcome model columns", {
  d <- .make_df_data(seed = 17L)

  expect_error(
    cate_dnn(
      data = d$df, x_cols = d$x_cols, y_col = "Y", t_col = "D",
      m_col = "m_oracle", m0_col = "m0_oracle", m1_col = "m1_oracle",
      e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "not both", ignore.case = TRUE
  )

  expect_error(
    cate_dnn(
      data = d$df, x_cols = d$x_cols, y_col = "Y", t_col = "D",
      e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "pseudo_outcome", ignore.case = TRUE
  )

  expect_error(
    cate_dnn(
      data = d$df, x_cols = d$x_cols, y_col = "Y", t_col = "D",
      m0_col = "m0_oracle", e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "pseudo_outcome", ignore.case = TRUE
  )

  expect_error(
    cate_dnn(
      data = d$df, x_cols = d$x_cols, y_col = "Y", t_col = "D",
      m1_col = "m1_oracle", e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "pseudo_outcome", ignore.case = TRUE
  )

  expect_error(
    cate_dnn(
      data = d$df, x_cols = d$x_cols, t_col = "D",
      m_col = "m_oracle", e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    )
  )
  expect_error(
    cate_dnn(
      data = d$df, x_cols = c("X1", "X_missing"),
      y_col = "Y", t_col = "D",
      m_col = "m_oracle", e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "specified by `x_cols`", ignore.case = TRUE
  )
})

test_that("data-frame API distinguishes NA vs non-finite treatment diagnostics", {
  d <- .make_df_data(seed = 19L)

  df_na <- d$df
  df_na$D[1] <- NA_real_
  expect_error(
    cate_dnn(
      data = df_na, x_cols = d$x_cols, y_col = "Y", t_col = "D",
      m_col = "m_oracle", e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "t_col.*NA",
    ignore.case = TRUE
  )

  df_inf <- d$df
  df_inf$D[1] <- Inf
  expect_error(
    cate_dnn(
      data = df_inf, x_cols = d$x_cols, y_col = "Y", t_col = "D",
      m_col = "m_oracle", e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "t_col.*non-finite",
    ignore.case = TRUE
  )
})

test_that("data-frame API distinguishes NA vs non-finite covariate diagnostics", {
  d <- .make_df_data(seed = 21L)
  base <- function(df) {
    tryCatch(
      cate_dnn(
        data = df, x_cols = d$x_cols, y_col = "Y", t_col = "D",
        m_col = "m_oracle", e_col = "e_oracle",
        X0 = d$X0, k = 2L, s = 3L, se = FALSE
      ),
      error = conditionMessage
    )
  }

  df_na <- d$df
  df_na[[d$x_cols[[1L]]]][2L] <- NA_real_
  df_inf <- d$df
  df_inf[[d$x_cols[[1L]]]][2L] <- Inf

  msg_na <- base(df_na)
  msg_inf <- base(df_inf)

  expect_match(msg_na, "X.*NA", ignore.case = TRUE)
  expect_match(msg_inf, "X.*non-finite", ignore.case = TRUE)
  expect_false(identical(msg_na, msg_inf))
})

test_that("data-frame API distinguishes NA vs non-finite x_cols index diagnostics", {
  d <- .make_df_data(seed = 29L)

  expect_error(
    cate_dnn(
      data = d$df, x_cols = c(1L, NA_integer_), y_col = "Y", t_col = "D",
      m_col = "m_oracle", e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "x_cols.*NA",
    ignore.case = TRUE
  )

  expect_error(
    cate_dnn(
      data = d$df, x_cols = c(1L, Inf), y_col = "Y", t_col = "D",
      m_col = "m_oracle", e_col = "e_oracle", X0 = d$X0, k = 2L, s = 3L
    ),
    regexp = "x_cols.*non-finite",
    ignore.case = TRUE
  )
})

# ============================================================================
# Interface-equivalence coverage: matrix vs data-frame paths
# ============================================================================

test_that("pseudo_outcome path: tau is identical across matrix and data interfaces", {
  d <- .make_df_data(seed = 101L)
  nms <- names(d$df)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])

  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw,
    k = 3L, s = 5L, method = "exact", se = FALSE
  )
  out_df_name <- cate_dnn(
    data = d$df, x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 3L, s = 5L, method = "exact", se = FALSE
  )
  out_df_idx <- cate_dnn(
    data = d$df, x_cols = match(d$x_cols, nms), pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 3L, s = 5L, method = "exact", se = FALSE
  )

  expect_equal(out_df_name$tau, out_mat$tau, tolerance = 1e-12)
  expect_equal(out_df_idx$tau, out_mat$tau, tolerance = 1e-12)
})

test_that("AIPW nuisance path matches matrix pseudo_outcome path", {
  d <- .make_df_data(seed = 102L)
  nms <- names(d$df)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])

  out_df_name <- cate_dnn(
    data = d$df,
    x_cols = d$x_cols,
    y_col = "Y",
    t_col = "D",
    m0_col = "m0_oracle",
    m1_col = "m1_oracle",
    e_col = "e_oracle",
    X0 = d$X0,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_df_idx <- cate_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    y_col = match("Y", nms),
    t_col = match("D", nms),
    m0_col = match("m0_oracle", nms),
    m1_col = match("m1_oracle", nms),
    e_col = match("e_oracle", nms),
    X0 = d$X0,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_df_name$tau, out_mat$tau, tolerance = 1e-12)
  expect_equal(out_df_idx$tau, out_mat$tau, tolerance = 1e-12)
})

test_that("Robinson nuisance path matches matrix pseudo_outcome path", {
  d <- .make_df_data(seed = 103L)
  nms <- names(d$df)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])

  out_df_name <- cate_dnn(
    data = d$df,
    x_cols = d$x_cols,
    y_col = "Y",
    t_col = "D",
    m_col = "m_oracle",
    e_col = "e_oracle",
    X0 = d$X0,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_df_idx <- cate_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    y_col = match("Y", nms),
    t_col = match("D", nms),
    m_col = match("m_oracle", nms),
    e_col = match("e_oracle", nms),
    X0 = d$X0,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_robinson,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_df_name$tau, out_mat$tau, tolerance = 1e-12)
  expect_equal(out_df_idx$tau, out_mat$tau, tolerance = 1e-12)
})

test_that("pseudo_outcome path: tau and se agree across interfaces (exact, fast jackknife)", {
  d <- .make_df_data(n = 45L, m = 4L, seed = 104L)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])

  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw,
    k = 2L, s = 4L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df <- cate_dnn(
    data = d$df, x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 2L, s = 4L, method = "exact", se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df$tau, out_mat$tau, tolerance = 1e-12)
  expect_equal(out_df$se, out_mat$se, tolerance = 1e-12)
})

test_that("pseudo_outcome path: tau and se agree across interfaces (exact, loo jackknife)", {
  d <- .make_df_data(n = 34L, m = 3L, seed = 105L)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])

  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw,
    k = 1L, s = 3L, method = "exact", se = TRUE, jackknife = "loo"
  )
  out_df <- cate_dnn(
    data = d$df, x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 1L, s = 3L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_equal(out_df$tau, out_mat$tau, tolerance = 1e-12)
  expect_equal(out_df$se, out_mat$se, tolerance = 1e-12)
})

test_that("pseudo_outcome path: tau and se agree across interfaces (asymp, fast jackknife)", {
  d <- .make_df_data(n = 48L, m = 3L, seed = 106L)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])

  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw,
    k = 2L, s = 5L, method = "asymp", se = TRUE, jackknife = "fast"
  )
  out_df <- cate_dnn(
    data = d$df, x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 2L, s = 5L, method = "asymp", se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df$tau, out_mat$tau, tolerance = 1e-12)
  expect_equal(out_df$se, out_mat$se, tolerance = 1e-12)
})

test_that("x0 path via data-frame pseudo_outcome matches matrix interface", {
  d <- .make_df_data(n = 40L, p = 3L, m = 5L, seed = 107L)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])
  j <- 3L

  out_mat <- cate_dnn(
    X = x_mat, x0 = d$X0[j, ], pseudo_outcome = d$psi_aipw,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_df <- cate_dnn(
    data = d$df, x_cols = d$x_cols, x0 = d$X0[j, ], pseudo_outcome = d$psi_aipw,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_df$tau, out_mat$tau, tolerance = 1e-12)
})

test_that("return extras are identical across interfaces on pseudo_outcome path", {
  d <- .make_df_data(n = 44L, m = 3L, seed = 108L)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])

  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw,
    k = 2L, s = 4L, method = "exact", se = FALSE,
    return = c("W", "ranks", "ord", "dx", "psi", "ess", "neighbors")
  )
  out_df <- cate_dnn(
    data = d$df, x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 2L, s = 4L, method = "exact", se = FALSE,
    return = c("W", "ranks", "ord", "dx", "psi", "ess", "neighbors")
  )

  expect_equal(out_df$W, out_mat$W, tolerance = 1e-12)
  expect_equal(out_df$ranks, out_mat$ranks)
  expect_equal(out_df$ord, out_mat$ord)
  expect_equal(out_df$dx, out_mat$dx, tolerance = 1e-12)
  expect_equal(out_df$psi, out_mat$psi, tolerance = 1e-12)
  expect_equal(out_df$ess, out_mat$ess, tolerance = 1e-12)
  expect_equal(out_df$neighbors, out_mat$neighbors)
})

test_that("meta fields are identical across interfaces on pseudo_outcome path", {
  d <- .make_df_data(seed = 109L)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])

  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_df <- cate_dnn(
    data = d$df, x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 2L, s = 4L, method = "exact", se = FALSE
  )

  expect_identical(out_df$meta$n, out_mat$meta$n)
  expect_identical(out_df$meta$m, out_mat$meta$m)
  expect_identical(out_df$meta$p, out_mat$meta$p)
  expect_identical(out_df$meta$k, out_mat$meta$k)
  expect_identical(out_df$meta$s, out_mat$meta$s)
  expect_identical(out_df$meta$method, out_mat$meta$method)
  expect_identical(out_df$meta$distance, out_mat$meta$distance)
})

test_that("distance spec object preserves interface equivalence", {
  d <- .make_df_data(seed = 110L)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])
  dist_spec <- distance_l2()

  out_mat <- cate_dnn(
    X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw,
    distance = dist_spec, k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_df <- cate_dnn(
    data = d$df, x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, distance = dist_spec, k = 2L, s = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_df$tau, out_mat$tau, tolerance = 1e-12)
})

test_that("tibble input matches data.frame input on pseudo_outcome path", {
  skip_if_not_installed("tibble")
  d <- .make_df_data(seed = 111L)

  out_df <- cate_dnn(
    data = d$df, x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_tbl <- cate_dnn(
    data = tibble::as_tibble(d$df), x_cols = d$x_cols, pseudo_outcome = d$psi_aipw,
    X0 = d$X0, k = 2L, s = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_tbl$tau, out_df$tau, tolerance = 1e-12)
})
