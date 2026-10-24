# Stored reporting-target contract for the shared plotting layer:
# every fitted result object stores the evaluation-point matrix `x0` so that
# downstream consumers (plotting, reporting) can reconstruct target profiles.

.x0_expectations <- function(fit, m, p, colnames_expected = NULL) {
  expect_true(is.matrix(fit$x0))
  expect_identical(dim(fit$x0), c(m, p))
  expect_true(is.numeric(fit$x0))
  if (!is.null(colnames_expected)) {
    expect_identical(colnames(fit$x0), colnames_expected)
  }
}

test_that("np_dnn and cate_dnn store x0 with propagated covariate names", {
  withr::local_preserve_seed()
  set.seed(21)
  n <- 60L
  df <- data.frame(
    X1 = rnorm(n), X2 = rnorm(n),
    Y = rnorm(n)
  )
  X0 <- rbind(c(0, 0), c(1, -1), c(-1, 1))

  fit_np <- np_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y",
    X0 = X0, k = 1L, s = 10L, se = FALSE
  )
  .x0_expectations(fit_np, 3L, 2L, colnames_expected = c("X1", "X2"))
  expect_equal(unname(fit_np$x0), X0)

  # Matrix interface without column names: x0 stored, no names invented.
  psi <- rnorm(n)
  X <- cbind(df$X1, df$X2)
  fit_cate <- cate_dnn(
    X = X, X0 = X0, pseudo_outcome = psi, k = 1L, s = 10L, se = FALSE
  )
  .x0_expectations(fit_cate, 3L, 2L)
  expect_null(colnames(fit_cate$x0))

  # Single-point x0 vector is stored as a 1 x p matrix.
  fit_one <- cate_dnn(
    X = X, x0 = c(0.5, -0.5), pseudo_outcome = psi, k = 1L, s = 10L,
    se = FALSE
  )
  .x0_expectations(fit_one, 1L, 2L)
  expect_equal(unname(fit_one$x0[1L, ]), c(0.5, -0.5))
})

test_that("user-supplied X0 column names are preserved over training names", {
  withr::local_preserve_seed()
  set.seed(22)
  n <- 50L
  df <- data.frame(X1 = rnorm(n), X2 = rnorm(n), Y = rnorm(n))
  X0 <- matrix(0, 2L, 2L, dimnames = list(NULL, c("a", "b")))

  fit <- np_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y",
    X0 = X0, k = 1L, s = 10L, se = FALSE
  )
  expect_identical(colnames(fit$x0), c("a", "b"))
})

test_that("plqr_dnn, gmm_dnn, and tdnn_np store x0", {
  withr::local_seed(123L)
  demo <- plqr_demo_data(n = 80L, p = 2L, p_d = 2L, tau = 0.5, seed = 5)
  df <- demo$data
  X <- as.matrix(df[paste0("X", 1:2)])
  D <- as.matrix(df[paste0("D", 1:2)])
  mD <- as.matrix(df[paste0("mD", 1:2)])
  X0 <- matrix(rnorm(4L), 2L, 2L)

  fit_plqr <- plqr_dnn(
    X = X, D = D, g = as.numeric(df$g_oracle), mD = mD,
    Y = as.numeric(df$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = X0, k = 1L, s = 8L, se = FALSE
  )
  .x0_expectations(fit_plqr, 2L, 2L, colnames_expected = c("X1", "X2"))

  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L, seed = 6L)
  model <- gmm_model_dr_iv(
    z_cols = paste0("Z", seq_len(dgp$q)),
    pi_cols = paste0("piZ", seq_len(dgp$q), "_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
  fit_gmm <- gmm_dnn(
    data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
    model = model, k = 1L, s = 12L
  )
  .x0_expectations(fit_gmm, 1L, 2L, colnames_expected = c("X1", "X2"))

  set.seed(23)
  n <- 70L
  Xm <- matrix(rnorm(n * 2L), n, 2L,
               dimnames = list(NULL, c("V1", "V2")))
  y <- rnorm(n)
  fit_tdnn <- np_tdnn(
    X = Xm, x0 = c(0, 0), Y = y, s1 = 4L, s2 = 8L, se = FALSE
  )
  .x0_expectations(fit_tdnn, 1L, 2L, colnames_expected = c("V1", "V2"))
})

test_that("stored x0 is invariant to se, jackknife, and return options", {
  withr::local_preserve_seed()
  set.seed(24)
  n <- 50L
  x_mat <- matrix(rnorm(n * 2L), n, 2L,
                  dimnames = list(NULL, c("X1", "X2")))
  y <- rnorm(n)
  X0 <- rbind(c(0, 0), c(1, 1))

  fit_a <- np_dnn(X = x_mat, X0 = X0, Y = y, k = 1L, s = 8L, se = TRUE)
  fit_b <- np_dnn(X = x_mat, X0 = X0, Y = y, k = 1L, s = 8L, se = FALSE,
                  return = c("W", "dx", "ess"))
  fit_c <- np_dnn(X = x_mat, X0 = X0, Y = y, k = 1L, s = 8L, se = TRUE,
                  jackknife = "loo")

  expect_identical(fit_a$x0, fit_b$x0)
  expect_identical(fit_a$x0, fit_c$x0)
  expect_identical(colnames(fit_a$x0), c("X1", "X2"))
})
