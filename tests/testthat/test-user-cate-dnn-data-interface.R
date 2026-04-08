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
