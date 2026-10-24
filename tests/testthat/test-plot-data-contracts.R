# Deeper contract and cross-consistency tests for the plotting data layer
# (Phase 0). Complements test-plot-data.R with:
#   - numerical identities against stored fit components,
#   - row-ordering contracts for multi-coefficient / multi-target frames,
#   - synthetic-object edge cases for stored-se shapes,
#   - exact zero-tail combinatorics of DNN weights,
#   - tie handling and determinism,
#   - signed TDNN weight profiles,
#   - cross-extractor consistency (weights frame vs diagnostics frame),
#   - CATE overlap and PLQR Gram diagnostic invariants.

.pdc_np_fit <- function(n = 60L, k = 1L, s = 10L, se = TRUE, return = c("W"),
                        X0 = rbind(c(0, 0), c(1, 0)), seed = 41L) {
  .test_rng_scope({
    x_mat <- matrix(rnorm(n * 2L), n, 2L,
                    dimnames = list(NULL, c("X1", "X2")))
    y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.3)
    np_dnn(
      X = x_mat, X0 = X0, Y = y_vec,
      k = k, s = s, se = se, return = return
    )
  }, seed = seed)
}

test_that("estimates frame reproduces stored components exactly", {
  .test_rng_scope({
    fit <- .pdc_np_fit(se = TRUE)
    df <- plot_data_estimates(fit)

    # Point estimates and ses are the stored values, in target order.
    expect_identical(df$estimate, as.numeric(fit$mu))
    expect_identical(df$se, as.numeric(fit$se))

    # Interval width scales exactly with the normal quantile ratio.
    df50 <- plot_data_estimates(fit, level = 0.5)
    ratio <- (df$ci_upper - df$ci_lower) / (df50$ci_upper - df50$ci_lower)
    expect_equal(ratio, rep(qnorm(0.975) / qnorm(0.75), nrow(df)))
    expect_identical(attr(df50, "level"), 0.5)

    # Extractors are pure: repeated calls give identical frames.
    expect_identical(plot_data_estimates(fit), df)

    # TDNN fits dispatch through the rdistnn_np_result method.
    set.seed(42)
    n <- 60L
    Xm <- matrix(rnorm(n * 2L), n, 2L)
    fit_tdnn <- np_tdnn(
      X = Xm, x0 = c(0, 0), Y = rnorm(n), s1 = 4L, s2 = 8L, se = FALSE
    )
    df_tdnn <- plot_data_estimates(fit_tdnn)
    expect_identical(unique(df_tdnn$coef), "mu")
    expect_identical(df_tdnn$estimate, as.numeric(fit_tdnn$mu))
  })
})

test_that("weight frames reject malformed retained distances before indexing", {
  fit <- .pdc_np_fit(n = 20L, s = 3L, se = FALSE, return = c("W", "dx"))
  for (value in list(fit$dx[-1L, , drop = FALSE], rep(1, length(fit$dx)),
      matrix(NA_real_, nrow(fit$dx), ncol(fit$dx)))) {
    bad <- fit
    bad$dx <- value
    err <- tryCatch(plot_data_weights(bad), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    if (inherits(err, "condition")) expect_match(conditionMessage(err), "Refit", fixed = TRUE)
  }
})

test_that("stored-se shapes are normalized and mismatches abort", {
  # Test the shape adapter directly; a fitted result must retain current
  # scientific state even when its optional reporting extras are omitted.
  base <- list(mu = c(0.1, 0.2, 0.3), se = NULL,
               x0 = matrix(seq_len(6L), 3L, 2L))
  extract <- function(fit) RdistNN:::.plot_data_estimates_frame(
    matrix(fit$mu, nrow = 1L), fit$se, fit$x0, "mu", .95)

  # Vector se and 1 x m matrix se give identical frames.
  fit_vec <- base
  fit_vec$se <- c(1, 2, 3)
  fit_mat <- base
  fit_mat$se <- matrix(c(1, 2, 3), nrow = 1L)
  expect_identical(
    extract(fit_vec),
    extract(fit_mat)
  )
  expect_identical(extract(fit_vec)$se, c(1, 2, 3))

  # Wrong dimensions abort with the dimension message.
  fit_bad <- base
  fit_bad$se <- matrix(1, nrow = 2L, ncol = 3L)
  expect_error(extract(fit_bad), "dimensions")

  # x0 without column names: coordinate columns, no name attribute.
  df <- extract(fit_vec)
  expect_true(all(c("x0_1", "x0_2") %in% names(df)))
  expect_null(attr(df, "x0_names"))
})

test_that("multi-coefficient frames interleave coefficients within targets", {
  demo <- plqr_demo_data(n = 80L, p = 2L, p_d = 2L, tau = 0.5, seed = 9)
  df_p <- demo$data
  fit <- plqr_dnn(
    X = as.matrix(df_p[paste0("X", 1:2)]),
    D = as.matrix(df_p[paste0("D", 1:2)]),
    g = as.numeric(df_p$g_oracle),
    mD = as.matrix(df_p[paste0("mD", 1:2)]),
    Y = as.numeric(df_p$Y),
    density = rep(2, nrow(df_p)),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = rbind(c(0, 0), c(0.5, -0.5)), k = 1L, s = 8L, se = TRUE
  )
  df <- plot_data_estimates(fit, quantity = "theta")
  est <- as.matrix(fit$theta)
  q <- nrow(est)

  # Row order: coefficients fastest, targets in blocks.
  expect_identical(df$target_id, rep(1:2, each = q))
  for (j in 1:2) {
    block <- df[df$target_id == j, ]
    expect_equal(block$estimate, as.numeric(est[, j]))
    expect_equal(block$se, as.numeric(as.matrix(fit$theta_se)[, j]))
    expect_identical(block$x0_1, rep(unname(fit$x0[j, 1L]), q))
  }
})

test_that("gmm frames map theta columns to targets with coordinates", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L,
                     seed = 10L)
  model <- gmm_model_dr_iv(
    z_cols = paste0("Z", seq_len(dgp$q)),
    pi_cols = paste0("piZ", seq_len(dgp$q), "_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
  X0 <- rbind(dgp$x0, dgp$x0 + 0.25)
  fit <- gmm_dnn(
    data = dgp$data, x_cols = c("X1", "X2"), X0 = X0,
    model = model, k = 1L, s = 12L
  )
  df <- plot_data_estimates(fit)
  p_theta <- length(fit$meta$coef_names)

  expect_identical(nrow(df), p_theta * 2L)
  expect_identical(df$coef, rep(fit$meta$coef_names, times = 2L))
  for (j in 1:2) {
    block <- df[df$target_id == j, ]
    expect_equal(block$estimate, as.numeric(as.matrix(fit$theta)[, j]))
    expect_equal(unique(block$x0_1), unname(X0[j, 1L]))
  }
})

test_that("weight frame values are a permutation with exact zero tail", {
  n <- 60L
  s <- 10L
  fit <- .pdc_np_fit(n = n, s = s, se = FALSE, return = c("W", "dx"))
  dw <- plot_data_weights(fit)

  for (j in seq_len(fit$meta$m)) {
    dj <- dw[dw$target_id == j, ]
    # Multiset identity with the stored weight column.
    expect_equal(sort(dj$weight), sort(fit$W[, j]))
    # Exact-weight combinatorics: support is n - s + 1 ranks, so the
    # distance-ordered tail contains exactly s - 1 zero weights.
    expect_identical(sum(dj$weight == 0), s - 1L)
    expect_true(all(dj$weight[seq_len(n - s + 1L)] > 0))
    # Cumulative share reaches one at the support edge already.
    expect_equal(dj$cum_share[n - s + 1L], 1)
    # First-rank weight is the DNN closed form s/n.
    expect_equal(dj$weight[1L], s / n)
  }
})

test_that("weight frame exact support uses n - s + k, not first s ranks", {
  n <- 60L
  k <- 3L
  s <- 10L
  fit <- .pdc_np_fit(n = n, k = k, s = s, se = FALSE, return = c("W", "dx"))
  dw <- plot_data_weights(fit)
  support_edge <- n - s + k

  for (j in seq_len(fit$meta$m)) {
    dj <- dw[dw$target_id == j, ]
    expect_identical(sum(dj$weight == 0), s - k)
    expect_true(all(dj$weight[seq_len(support_edge)] > 0))
    expect_true(all(dj$weight[(support_edge + 1L):n] == 0))
    expect_equal(dj$cum_share[support_edge], 1)
    expect_false(identical(support_edge, s))
  }
})

test_that("distance ties are ordered deterministically", {
  .test_rng_scope({
    set.seed(43)
    n <- 40L
    x_base <- matrix(rnorm(n * 2L), n, 2L)
    # Force exact distance ties: duplicate a block of rows.
    x_base[21:30, ] <- x_base[1:10, ]
    y <- rnorm(n)
    expect_warning(fit <- np_dnn(
      X = x_base, x0 = c(0, 0), Y = y,
      k = 1L, s = 6L, se = FALSE, return = c("W", "dx")
    ), "Exact distance ties detected", fixed = TRUE)
    expect_warning(dw1 <- plot_data_weights(fit),
                   "Exact distance ties detected", fixed = TRUE)
    expect_warning(dw2 <- plot_data_weights(fit),
                   "Exact distance ties detected", fixed = TRUE)

    expect_identical(dw1, dw2)
    expect_identical(sort(dw1$obs_id), seq_len(n))
    expect_true(!is.unsorted(dw1$distance))
    # Within exact distance ties, ranks are ordered by descending weight.
    dup <- duplicated(dw1$distance) | duplicated(dw1$distance, fromLast = TRUE)
    if (any(dup)) {
      for (d in unique(dw1$distance[dup])) {
        w_tie <- dw1$weight[dw1$distance == d]
        expect_true(!is.unsorted(rev(w_tie)))
      }
    }
  })
})

test_that("TDNN weight frames handle signed weights", {
  .test_rng_scope({
    set.seed(44)
    n <- 60L
    Xm <- matrix(rnorm(n * 2L), n, 2L)
    fit <- np_tdnn(
      X = Xm, x0 = c(0, 0), Y = rnorm(n),
      s1 = 4L, s2 = 8L, se = FALSE, return = c("W", "dx")
    )
    dw <- plot_data_weights(fit)

    # Two-scale weights are signed but still normalized.
    expect_true(any(dw$weight < 0))
    expect_equal(sum(dw$weight), 1)
    expect_equal(dw$cum_share[nrow(dw)], 1)
    expect_true(!is.unsorted(dw$distance))
    expect_identical(sort(dw$obs_id), seq_len(n))
  })
})

test_that("weights frame and diagnostics frame agree target by target", {
  fit <- .pdc_np_fit(se = FALSE, return = c("W", "dx"))
  dw <- plot_data_weights(fit)
  dd <- plot_data_diagnostics(fit)

  for (j in seq_len(fit$meta$m)) {
    dj <- dw[dw$target_id == j, ]
    expect_equal(dd$n_eff[j], 1 / sum(dj$weight^2))
    expect_equal(dd$max_weight[j], max(dj$weight))
    expect_equal(dd$top1_weight_share[j], max(dj$weight) / sum(dj$weight))
    # Every represented nonzero weight belongs to diagnostic support,
    # including tail weights below sqrt(.Machine$double.eps).
    expect_equal(
      dd$n_nonzero_weight[j],
      sum(dj$weight != 0)
    )
    # The rank-s radius equals the distance at rank s in the weight frame.
    expect_equal(dd$radius_rank_s[j], dj$distance[fit$meta$s])
  }
})

test_that("CATE overlap diagnostics satisfy mass and propensity invariants", {
  .test_rng_scope({
    set.seed(45)
    n <- 120L
    df <- data.frame(
      X1 = rnorm(n), X2 = rnorm(n)
    )
    e_true <- plogis(0.5 * df$X1)
    df$T <- as.integer(runif(n) < e_true)
    df$e <- e_true
    df$m0 <- df$X1
    df$m1 <- df$X1 + 1
    df$Y <- ifelse(df$T == 1L, df$m1, df$m0) + rnorm(n, sd = 0.3)

    fit <- cate_dnn(
      data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
      e_col = "e", m0_col = "m0", m1_col = "m1",
      X0 = rbind(c(0, 0), c(1, 0)), k = 1L, s = 15L, se = FALSE
    )
    dd <- plot_data_diagnostics(fit)

    expect_true(all(dd$overlap_available))
    # Simplex weights split exactly into treated and control mass.
    expect_equal(dd$treated_mass + dd$control_mass, rep(1, nrow(dd)))
    expect_true(all(dd$treated_share >= 0 & dd$treated_share <= 1))
    expect_equal(dd$treated_share, dd$treated_mass)
    # Propensity summaries are ordered and interior.
    expect_true(all(dd$propensity_available))
    expect_true(all(dd$weighted_propensity_q10 <= dd$weighted_mean_propensity))
    expect_true(all(dd$weighted_mean_propensity <= dd$weighted_propensity_q90))
    expect_true(all(dd$weighted_mean_propensity > 0 &
                      dd$weighted_mean_propensity < 1))

    # Matrix-interface fit stores no treatment: overlap unavailable, frame
    # still well-formed.
    fit_mat <- cate_dnn(
      X = as.matrix(df[, c("X1", "X2")]),
      X0 = matrix(c(0, 0), nrow = 1L),
      pseudo_outcome = rnorm(n), k = 1L, s = 15L, se = FALSE
    )
    dd_mat <- plot_data_diagnostics(fit_mat)
    expect_false(any(dd_mat$overlap_available))
    expect_true(all(is.na(dd_mat$treated_mass)))
  })
})

test_that("PLQR Gram diagnostics satisfy eigenvalue invariants", {
  demo <- plqr_demo_data(n = 90L, p = 2L, p_d = 2L, tau = 0.5, seed = 11)
  df_p <- demo$data
  fit <- plqr_dnn(
    X = as.matrix(df_p[paste0("X", 1:2)]),
    D = as.matrix(df_p[paste0("D", 1:2)]),
    g = as.numeric(df_p$g_oracle),
    mD = as.matrix(df_p[paste0("mD", 1:2)]),
    Y = as.numeric(df_p$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = rbind(c(0, 0), c(0.5, 0.5)), k = 1L, s = 8L, se = FALSE
  )
  dd <- plot_data_diagnostics(fit)

  expect_true(all(c("gram_available", "gram_min_eig",
                    "gram_condition_number") %in% names(dd)))
  avail <- which(dd$gram_available)
  expect_true(length(avail) > 0)
  expect_true(all(dd$gram_min_eig[avail] <= dd$gram_max_eig[avail]))
  expect_true(all(dd$gram_min_eig[avail] >= 0))
  expect_equal(
    dd$gram_condition_number[avail],
    dd$gram_max_eig[avail] / dd$gram_min_eig[avail]
  )
})
