# Tests for the plotting data layer (Phase 0):
# plot_data_estimates(), plot_data_weights(), plot_data_diagnostics().
# See markdown_notes/design/13_estimator_and_result_architecture.md.

.pd_np_fit <- function(n = 60L, se = TRUE, return = c("W"), X0 = NULL) {
  set.seed(31)
  x_mat <- matrix(rnorm(n * 2L), n, 2L,
                  dimnames = list(NULL, c("X1", "X2")))
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.3)
  if (is.null(X0)) X0 <- rbind(c(0, 0), c(1, 0), c(-1, 0))
  np_dnn(
    X = x_mat, X0 = X0, Y = y_vec,
    k = 1L, s = 10L, se = se, return = return
  )
}

test_that("plot_data_estimates returns the documented frame for np fits", {
  fit <- .pd_np_fit(se = TRUE)
  df <- plot_data_estimates(fit)

  expect_s3_class(df, "data.frame")
  expect_identical(
    names(df),
    c("target_id", "coef", "estimate", "se", "ci_lower", "ci_upper",
      "x0_1", "x0_2")
  )
  expect_identical(df$target_id, 1:3)
  expect_identical(unique(df$coef), "mu")
  expect_equal(df$estimate, as.numeric(fit$mu))
  z <- qnorm(0.975)
  expect_equal(df$ci_lower, df$estimate - z * df$se)
  expect_equal(df$ci_upper, df$estimate + z * df$se)
  expect_identical(attr(df, "x0_names"), c("X1", "X2"))
  expect_identical(attr(df, "level"), 0.95)
  expect_equal(df$x0_1, c(0, 1, -1))

  # Narrower level gives narrower bounds.
  df80 <- plot_data_estimates(fit, level = 0.8)
  expect_true(all(df80$ci_upper - df80$ci_lower <
                    df$ci_upper - df$ci_lower))
  expect_error(plot_data_estimates(fit, level = 1.2), "in \\(0,1\\)")
})

test_that("plot_data_estimates keeps disabled uncertainty and rejects incomplete targets", {
  fit <- .pd_np_fit(se = FALSE)
  df <- plot_data_estimates(fit)
  expect_true(all(is.na(df$se)))
  expect_true(all(is.na(df$ci_lower)))

  # Current snapshots retain their targets even without uncertainty.
  fit$x0 <- NULL
  expect_error(plot_data_estimates(fit), "x0.*[Rr]efit")
})

test_that("plot_data_estimates dispatches for cate, plqr, and gmm fits", {
  set.seed(32)
  n <- 60L
  X <- matrix(rnorm(n * 2L), n, 2L)
  psi <- rnorm(n, mean = X[, 1L])
  fit_cate <- cate_dnn(
    X = X, X0 = rbind(c(0, 0), c(1, 1)), pseudo_outcome = psi,
    k = 1L, s = 10L, se = FALSE
  )
  df_cate <- plot_data_estimates(fit_cate)
  expect_identical(unique(df_cate$coef), "tau")
  expect_identical(nrow(df_cate), 2L)

  demo <- plqr_demo_data(n = 80L, p = 2L, p_d = 2L, tau = 0.5, seed = 7)
  df_p <- demo$data
  fit_plqr <- plqr_dnn(
    X = as.matrix(df_p[paste0("X", 1:2)]),
    D = as.matrix(df_p[paste0("D", 1:2)]),
    g = as.numeric(df_p$g_oracle),
    mD = as.matrix(df_p[paste0("mD", 1:2)]),
    Y = as.numeric(df_p$Y),
    density = rep(2, nrow(df_p)),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = matrix(0, 2L, 2L), k = 1L, s = 8L, se = TRUE
  )
  df_plqr <- plot_data_estimates(fit_plqr)
  p_d <- nrow(as.matrix(fit_plqr$theta))
  expect_identical(nrow(df_plqr), p_d * 2L)
  expect_identical(length(unique(df_plqr$coef)), p_d)
  expect_identical(df_plqr$target_id, rep(1:2, each = p_d))
  # Default PLQR quantity is the localized moment (no root solver exists);
  # this analytic-meat fit stores theta_se but no moment_se.
  expect_identical(attr(df_plqr, "quantity"), "moment")
  expect_equal(df_plqr$estimate, as.numeric(as.matrix(fit_plqr$moment)))
  expect_true(all(is.na(df_plqr$se)))

  df_plqr_theta <- plot_data_estimates(fit_plqr, quantity = "theta")
  expect_identical(attr(df_plqr_theta, "quantity"), "theta")
  expect_equal(df_plqr_theta$estimate, as.numeric(as.matrix(fit_plqr$theta)))
  expect_false(anyNA(df_plqr_theta$se))

  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L,
                     seed = 8L)
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
  df_gmm <- plot_data_estimates(fit_gmm)
  expect_identical(unique(df_gmm$coef), fit_gmm$meta$coef_names)
  expect_identical(nrow(df_gmm), length(fit_gmm$meta$coef_names))

  expect_error(plot_data_estimates(lm(y ~ x, data = data.frame(y = 1:5, x = 1:5))),
               "not implemented")
})

test_that("plot_data_weights orders ranks and computes cumulative shares", {
  fit <- .pd_np_fit(se = FALSE, return = c("W", "dx"))
  dw <- plot_data_weights(fit)

  n <- fit$meta$n
  m <- fit$meta$m
  expect_identical(
    names(dw),
    c("target_id", "rank", "obs_id", "weight", "cum_share", "distance")
  )
  expect_identical(nrow(dw), n * m)

  for (j in seq_len(m)) {
    dj <- dw[dw$target_id == j, ]
    expect_identical(dj$rank, seq_len(n))
    expect_identical(sort(dj$obs_id), seq_len(n))
    # Distance-based ranks: distances nondecreasing in rank.
    expect_true(!is.unsorted(dj$distance))
    # Normalized weights: cumulative share reaches one.
    expect_equal(dj$cum_share[n], 1)
    expect_equal(sum(dj$weight), 1)
  }

  # Without stored distances: ranks by descending weight, distance is NA.
  fit_nodx <- .pd_np_fit(se = FALSE, return = c("W"))
  dw2 <- plot_data_weights(fit_nodx)
  expect_true(all(is.na(dw2$distance)))
  d1 <- dw2[dw2$target_id == 1L, ]
  expect_true(!is.unsorted(rev(d1$weight)))

  # Missing W: actionable error naming the remedy.
  fit_now <- .pd_np_fit(se = FALSE, return = c("ess"))
  expect_error(plot_data_weights(fit_now), class = "rdistnn_extra_unavailable")
})

test_that("plot_data_diagnostics tidies the diagnostics frame", {
  fit <- .pd_np_fit(se = FALSE, return = c("W", "dx"))
  dd <- plot_data_diagnostics(fit)

  expect_identical(class(dd), "data.frame")
  expect_identical(dd$target_id, seq_len(nrow(dd)))
  expect_true(all(c("target_label", "n_eff", "max_weight",
                    "top5_weight_share") %in% names(dd)))
  expect_identical(dd$target_label, paste0("target_", seq_len(nrow(dd))))
  expect_true(all(dd$n_eff > 0))
})
