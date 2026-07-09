# Guardrail tests for the default kernel order `s`.
#
# History: the user-level wrappers originally shipped with a silent `s = 1L`
# default, which yields exactly uniform DNN weights and therefore no
# localization at all (every target point receives the global weighted mean).
# These tests pin the corrected behavior: when both `s` and `tuning` are
# omitted, `s` must resolve to the paper undersmoothing rule
# s = floor(n^((kx + 1)/(kx + 2))) via tuning_dnn_paper(), matching the
# default implemented in the Stata companion (cate_dnn.ado).  They are the
# `s`-analogue of test-k1-documentation-defaults.R.

.paper_s <- function(n, kx) as.integer(floor(n^((kx + 1) / (kx + 2))))

test_that("np_dnn default s follows the paper rule and records provenance", {
  set.seed(11)
  n <- 120L
  kx <- 2L
  X <- matrix(rnorm(n * kx), n, kx)
  y <- rnorm(n, mean = X[, 1L])

  fit <- np_dnn(X = X, x0 = c(0, 0), Y = y, se = FALSE, return = "tuning")

  expect_identical(fit$meta$s, .paper_s(n, kx))
  expect_identical(fit$meta$s, tuning_eval("paper_undersmooth", n = n, kx = kx, k = 1L)$s)
  expect_identical(fit$meta$tuning$source, "paper_default")
})

test_that("cate_dnn default s follows the paper rule and is localized", {
  set.seed(12)
  n <- 150L
  kx <- 2L
  X <- matrix(rnorm(n * kx), n, kx)
  psi <- rnorm(n, mean = 2 * X[, 1L])
  X0 <- rbind(c(-2, 0), c(2, 0))

  fit <- cate_dnn(X = X, X0 = X0, pseudo_outcome = psi, se = FALSE,
                  return = "tuning")

  expect_identical(fit$meta$s, .paper_s(n, kx))
  expect_identical(fit$meta$tuning$source, "paper_default")
  # Localization sanity: estimates at x1 = -2 and x1 = +2 must differ.
  # (Under the old s = 1 default both equalled the global mean exactly.)
  expect_gt(abs(fit$tau[2L] - fit$tau[1L]), 1)
})

test_that("plqr_dnn and gmm_dnn default s follows the paper rule", {
  demo <- plqr_demo_data(n = 90L, p = 2L, p_d = 2L, tau = 0.5, seed = 3)
  df <- demo$data
  X <- as.matrix(df[paste0("X", 1:2)])
  D <- as.matrix(df[paste0("D", 1:2)])
  mD <- as.matrix(df[paste0("mD", 1:2)])
  fit_plqr <- plqr_dnn(
    X = X, D = D, g = as.numeric(df$g_oracle), mD = mD,
    Y = as.numeric(df$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = matrix(0, 1L, 2L), se = FALSE, return = "tuning"
  )
  expect_identical(fit_plqr$meta$s, .paper_s(nrow(X), ncol(X)))
  expect_identical(fit_plqr$meta$tuning$source, "paper_default")

  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 120L, kx = 2L, q = 2L, seed = 4L)
  model <- gmm_model_dr_iv(
    z_cols = paste0("Z", seq_len(dgp$q)),
    pi_cols = paste0("piZ", seq_len(dgp$q), "_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
  fit_gmm <- gmm_dnn(
    data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
    model = model, return = "tuning"
  )
  expect_identical(fit_gmm$meta$s, .paper_s(120L, 2L))
  expect_identical(fit_gmm$meta$tuning$source, "paper_default")
})

test_that("explicit s is respected and marked manual", {
  set.seed(13)
  n <- 80L
  X <- matrix(rnorm(n * 2L), n, 2L)
  y <- rnorm(n)
  fit <- np_dnn(X = X, x0 = c(0, 0), Y = y, s = 7L, se = FALSE,
                return = "tuning")
  expect_identical(fit$meta$s, 7L)
  # Manual s does not attach a tuning payload to meta (tuning_used = FALSE).
  expect_null(fit$meta$tuning)
})

test_that("degenerate (near-uniform) weights trigger a warning", {
  set.seed(14)
  n <- 60L
  X <- matrix(rnorm(n * 2L), n, 2L)
  y <- rnorm(n)

  # s = 1 gives exactly uniform weights.
  expect_warning(
    np_dnn(X = X, x0 = c(0, 0), Y = y, k = 1L, s = 1L, se = FALSE),
    "near-.?uniform|localization"
  )
  # k = s also gives exactly uniform weights.
  expect_warning(
    np_dnn(X = X, x0 = c(0, 0), Y = y, k = 3L, s = 3L, se = FALSE),
    "near-.?uniform|localization"
  )
  # A healthy configuration must stay silent.
  expect_no_warning(
    np_dnn(X = X, x0 = c(0, 0), Y = y, k = 1L, s = 15L, se = FALSE)
  )
})

test_that("R default agrees with the Stata companion default formula", {
  # RdistNN_stata/src/cate_dnn.ado uses: floor(n^((kx + 1)/(kx + 2))).
  for (n in c(100L, 250L, 500L, 1000L)) {
    for (kx in c(1L, 2L, 3L)) {
      expect_identical(
        tuning_eval("paper_undersmooth", n = n, kx = kx, k = 1L)$s,
        as.integer(floor(n^((kx + 1) / (kx + 2))))
      )
    }
  }
})
