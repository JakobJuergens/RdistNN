skip_if_not(
  all(vapply(
    c("np_dnn", "diagnostics"),
    exists, logical(1),
    mode = "function"
  )),
  "np_dnn() or diagnostics() not found; skipping file"
)

set.seed(777)

.diag_top_share <- function(w, k) {
  k_use <- max(1L, min(length(w), as.integer(k)))
  sum(sort(as.numeric(w), decreasing = TRUE)[seq_len(k_use)])
}

.diag_weighted_median <- function(x, w) {
  .diag_weighted_quantile(x, w, prob = 0.5)
}

.diag_weighted_quantile <- function(x, w, prob = 0.5) {
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- w[ord]
  cw <- cumsum(w_ord) / sum(w_ord)
  x_ord[which(cw >= prob)[1L]]
}

.diag_weighted_sd <- function(x, w) {
  mu <- sum(x * w) / sum(w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}

test_that("diagnostics(np_result) returns expected schema and identities", {
  n <- 45L
  p <- 2L
  m <- 3L
  x_mat <- matrix(rnorm(n * p), n, p)
  x0_mat <- matrix(rnorm(m * p), m, p)
  y_vec <- x_mat[, 1L] - 0.4 * x_mat[, 2L] + rnorm(n, sd = 0.2)

  fit <- np_dnn(
    X = x_mat, X0 = x0_mat, Y = y_vec,
    distance = "l2", method = "exact",
    k = 2L, s = 4L, se = FALSE,
    return = c("W", "dx")
  )
  d <- diagnostics(fit)

  expect_s3_class(d, "dnn_diagnostics")
  expect_equal(nrow(d), m)
  expect_true(all(c(
    "target_id", "n_eff", "max_weight",
    "top1_weight_share", "top5_weight_share", "top10_weight_share",
    "radius_rank_k", "radius_rank_s",
    "weighted_mean_distance", "weighted_median_distance",
    "outcome_available", "weighted_mean_y", "weighted_sd_y",
    "weighted_y_q10", "weighted_y_q50", "weighted_y_q90"
  ) %in% names(d)))

  expect_equal(d$n_eff, 1 / colSums(fit$W^2), tolerance = 1e-12)
  expect_equal(d$max_weight, apply(fit$W, 2L, max), tolerance = 1e-12)
  expect_equal(
    d$top1_weight_share,
    apply(fit$W, 2L, .diag_top_share, k = 1L),
    tolerance = 1e-12
  )
  expect_equal(
    d$top5_weight_share,
    apply(fit$W, 2L, .diag_top_share, k = 5L),
    tolerance = 1e-12
  )
  expect_equal(
    d$top10_weight_share,
    apply(fit$W, 2L, .diag_top_share, k = 10L),
    tolerance = 1e-12
  )

  k_use <- fit$meta$k
  s_use <- fit$meta$s
  expect_equal(
    d$radius_rank_k,
    apply(fit$dx, 2L, function(col) sort(col, partial = k_use)[k_use]),
    tolerance = 1e-12
  )
  expect_equal(
    d$radius_rank_s,
    apply(fit$dx, 2L, function(col) sort(col, partial = s_use)[s_use]),
    tolerance = 1e-12
  )
  expect_equal(
    d$weighted_mean_distance,
    colSums(fit$dx * fit$W) / colSums(fit$W),
    tolerance = 1e-12
  )
  expect_equal(
    d$weighted_median_distance,
    vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_median(fit$dx[, j], fit$W[, j]), numeric(1L)),
    tolerance = 1e-12
  )
  expect_true(all(d$outcome_available))
  expect_equal(d$weighted_mean_y, fit$mu, tolerance = 1e-12)
  expect_equal(
    d$weighted_sd_y,
    vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_sd(y_vec, fit$W[, j]), numeric(1L)),
    tolerance = 1e-12
  )
  expect_equal(
    d$weighted_y_q10,
    vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_quantile(y_vec, fit$W[, j], prob = 0.1), numeric(1L)),
    tolerance = 1e-12
  )
  expect_equal(
    d$weighted_y_q50,
    vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_median(y_vec, fit$W[, j]), numeric(1L)),
    tolerance = 1e-12
  )
  expect_equal(
    d$weighted_y_q90,
    vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_quantile(y_vec, fit$W[, j], prob = 0.9), numeric(1L)),
    tolerance = 1e-12
  )
})

test_that("general diagnostics reject unsupported objects and malformed retained matrices", {
  n <- 35L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)
  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = c("W", "dx")
  )

  expect_error(diagnostics(list()), "not implemented for class 'list'")

  bad_w <- fit
  bad_w$W[1L, 1L] <- NA_real_
  expect_error(diagnostics(bad_w), "W cannot contain `NA` values")

  bad_dx_dim <- fit
  bad_dx_dim$dx <- fit$dx[-1L, , drop = FALSE]
  expect_error(diagnostics(bad_dx_dim), "same dimensions as `W`")

  bad_dx_value <- fit
  bad_dx_value$dx[1L, 1L] <- Inf
  expect_error(diagnostics(bad_dx_value), "dx cannot contain non-finite values")
})

test_that("diagnostics(np_result) reports unavailable distance fields when dx is absent", {
  n <- 40L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)

  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = "W"
  )
  d <- diagnostics(fit)

  expect_true(all(!d$distance_diagnostics_available))
  expect_true(all(is.na(d$radius_rank_k)))
  expect_true(all(is.na(d$radius_rank_s)))
  expect_true(all(is.na(d$weighted_mean_distance)))
  expect_true(all(is.na(d$weighted_median_distance)))
})

test_that("diagnostics(np_result) errors clearly when W is not retained", {
  n <- 35L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)

  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = "ess"
  )

  expect_error(
    diagnostics(fit),
    "refit with `return` including \"W\""
  )
})

test_that("diagnostics(np_result) marks outcome unavailable for legacy objects without outcome state", {
  n <- 35L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)

  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = "W"
  )
  fit$.model <- NULL
  d <- diagnostics(fit)

  expect_true(all(!d$outcome_available))
  expect_true(all(is.na(d$weighted_mean_y)))
  expect_true(all(is.na(d$weighted_sd_y)))
})

test_that("diagnostics(np_result) validates retained outcome state", {
  n <- 35L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)
  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = "W"
  )

  bad_length <- fit
  bad_length$.model$y_vec <- y_vec[-1L]
  expect_error(diagnostics(bad_length), "length equal to nrow\\(object\\$W\\)")

  bad_na <- fit
  bad_na$.model$y_vec[1L] <- NA_real_
  expect_error(diagnostics(bad_na), "object\\$\\.model\\$y_vec cannot contain `NA` values")
})

test_that("diagnostics(np_result) includes tuning metadata summary", {
  demo <- np_demo_data(n = 80L, p = 2L, seed = 19L)
  x_mat <- as.matrix(demo$data[, c("X1", "X2")])
  y_vec <- demo$data$Y
  tun <- tuning_dnn_paper(kappa = 1.25)

  fit <- np_dnn(
    X = x_mat, x0 = demo$x0, Y = y_vec,
    method = "exact", k = 2L, tuning = tun, se = FALSE,
    return = "W"
  )
  d <- diagnostics(fit)

  expect_true(all(d$tuning_source == "tuning"))
  expect_true(all(d$tuning_name == "paper_undersmooth"))
  expect_true(all(d$tuning_grid_available))
  expect_true(all(d$tuning_default_s == fit$meta$s))
  expect_true(all(d$tuning_grid_n_candidates == nrow(attr(d, "tuning_grid"))))
  expect_true(all(d$tuning_kappa_default == 1.25))
})

test_that("diagnostics(cate_result) reports overlap when treatment state is available", {
  demo <- cate_demo_data(n = 70L, p = 2L, design = "hetero", sigma = 0.2, seed = 20L)
  fit <- cate_dnn(
    data = demo$data,
    x_cols = c("X1", "X2"),
    y_col = "Y",
    t_col = "D",
    m_col = "m_oracle",
    e_col = "e_oracle",
    x0 = demo$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )
  d <- diagnostics(fit)

  expect_s3_class(d, "dnn_diagnostics")
  expect_true(all(d$overlap_available))
  expect_true(all(d$propensity_available))

  t_vec <- demo$data$D
  e <- demo$data$e_oracle
  expect_equal(d$treated_mass, as.numeric(crossprod(t_vec, fit$W)), tolerance = 1e-12)
  expect_equal(d$control_mass, as.numeric(crossprod(1 - t_vec, fit$W)), tolerance = 1e-12)
  expect_equal(
    d$weighted_mean_propensity,
    as.numeric(crossprod(e, fit$W)) / colSums(fit$W),
    tolerance = 1e-12
  )
  expect_equal(d$treated_share, d$treated_mass / (d$treated_mass + d$control_mass), tolerance = 1e-12)
  expect_equal(
    d$weighted_propensity_q10,
    vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_quantile(e, fit$W[, j], prob = 0.1), numeric(1L)),
    tolerance = 1e-12
  )
  expect_equal(
    d$weighted_propensity_q90,
    vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_quantile(e, fit$W[, j], prob = 0.9), numeric(1L)),
    tolerance = 1e-12
  )
})

test_that("diagnostics(cate_result) separates overlap and propensity availability", {
  demo <- cate_demo_data(n = 60L, p = 2L, design = "hetero", sigma = 0.2, seed = 22L)
  fit <- cate_dnn(
    data = demo$data,
    x_cols = c("X1", "X2"),
    y_col = "Y",
    t_col = "D",
    m_col = "m_oracle",
    e_col = "e_oracle",
    x0 = demo$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )
  fit$.model$e <- NULL
  d <- diagnostics(fit)

  expect_true(all(d$overlap_available))
  expect_true(all(!d$propensity_available))
  expect_false(any(is.na(d$treated_mass)))
  expect_true(all(is.na(d$weighted_mean_propensity)))
  expect_true(all(is.na(d$weighted_propensity_q10)))
  expect_true(all(is.na(d$weighted_propensity_q90)))
})

test_that("diagnostics(cate_result) validates retained treatment and propensity state", {
  demo <- cate_demo_data(n = 60L, p = 2L, design = "hetero", sigma = 0.2, seed = 23L)
  fit <- cate_dnn(
    data = demo$data,
    x_cols = c("X1", "X2"),
    y_col = "Y",
    t_col = "D",
    m_col = "m_oracle",
    e_col = "e_oracle",
    x0 = demo$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )

  bad_length <- fit
  bad_length$.model$t_vec <- bad_length$.model$t_vec[-1L]
  expect_error(diagnostics(bad_length), "length equal to nrow\\(object\\$W\\)")

  bad_binary <- fit
  bad_binary$.model$t_vec[1L] <- 2
  expect_error(diagnostics(bad_binary), "binary in \\{0,1\\}")

  bad_e_length <- fit
  bad_e_length$.model$e <- bad_e_length$.model$e[-1L]
  expect_error(diagnostics(bad_e_length), "object\\$\\.model\\$e.*length equal")

  bad_e_na <- fit
  bad_e_na$.model$e[1L] <- NA_real_
  expect_error(diagnostics(bad_e_na), "object\\$\\.model\\$e cannot contain `NA` values")
})

test_that("diagnostics(cate_result) marks overlap unavailable on pseudo-outcome-only path", {
  n <- 45L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  fit <- cate_dnn(
    X = x_mat,
    x0 = c(0, 0),
    pseudo_outcome = rnorm(n),
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )
  d <- diagnostics(fit)

  expect_true(all(!d$overlap_available))
  expect_true(all(!d$propensity_available))
  expect_true(all(is.na(d$treated_mass)))
  expect_true(all(is.na(d$weighted_mean_propensity)))
})

test_that("diagnostics(plqr_result) reports residual-regressor Gram diagnostics", {
  demo <- plqr_demo_data(n = 75L, p = 2L, p_d = 1L, tau = 0.5, seed = 21L)
  fit <- plqr_dnn(
    data = demo$data,
    x_cols = c("X1", "X2"),
    D_cols = "D1",
    q_col = "q_oracle",
    mD_cols = "mD1",
    y_col = "Y",
    beta = demo$truth$beta_true,
    tau = 0.5,
    x0 = demo$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )
  d <- diagnostics(fit)

  d_resid <- demo$data$D1 - demo$data$mD1
  expected_gram <- as.numeric(crossprod(d_resid^2, fit$W))

  expect_s3_class(d, "dnn_diagnostics")
  expect_true(all(d$gram_available))
  expect_equal(d$gram_min_eig, expected_gram, tolerance = 1e-12)
  expect_equal(d$gram_max_eig, expected_gram, tolerance = 1e-12)
  expect_equal(d$gram_condition_number, rep(1, nrow(d)), tolerance = 1e-12)

  old_fit <- fit
  old_fit$.model <- NULL
  d_old <- diagnostics(old_fit)
  expect_true(all(!d_old$gram_available))
  expect_true(all(is.na(d_old$gram_min_eig)))
})

test_that("diagnostics(plqr_result) reports multi-column and singular Gram diagnostics", {
  demo <- plqr_demo_data(n = 85L, p = 2L, p_d = 2L, tau = 0.5, seed = 24L)
  fit <- plqr_dnn(
    data = demo$data,
    x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"),
    q_col = "q_oracle",
    mD_cols = c("mD1", "mD2"),
    y_col = "Y",
    beta = demo$truth$beta_true,
    tau = 0.5,
    x0 = demo$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )
  d <- diagnostics(fit)

  d_resid <- fit$.model$d_resid_mat
  expected <- lapply(seq_len(ncol(fit$W)), function(j) {
    G <- crossprod(d_resid, d_resid * as.numeric(fit$W[, j]))
    G <- (G + t(G)) / 2
    eigen(G, symmetric = TRUE, only.values = TRUE)$values
  })
  expect_equal(d$gram_min_eig, vapply(expected, min, numeric(1L)), tolerance = 1e-12)
  expect_equal(d$gram_max_eig, vapply(expected, max, numeric(1L)), tolerance = 1e-12)
  expect_true(all(is.finite(d$gram_condition_number)))

  singular <- fit
  singular$.model$d_resid_mat[, 2L] <- singular$.model$d_resid_mat[, 1L]
  d_singular <- diagnostics(singular)
  expect_true(all(d_singular$gram_min_eig <= d_singular$gram_rank_tol))
  expect_true(all(is.infinite(d_singular$gram_condition_number)))
})

test_that("diagnostics(plqr_result) validates retained residual-regressor state", {
  demo <- plqr_demo_data(n = 70L, p = 2L, p_d = 2L, tau = 0.5, seed = 25L)
  fit <- plqr_dnn(
    data = demo$data,
    x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"),
    q_col = "q_oracle",
    mD_cols = c("mD1", "mD2"),
    y_col = "Y",
    beta = demo$truth$beta_true,
    tau = 0.5,
    x0 = demo$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )

  bad_rows <- fit
  bad_rows$.model$d_resid_mat <- bad_rows$.model$d_resid_mat[-1L, , drop = FALSE]
  expect_error(diagnostics(bad_rows), "row count equal to nrow\\(object\\$W\\)")

  bad_value <- fit
  bad_value$.model$d_resid_mat[1L, 1L] <- NaN
  expect_error(diagnostics(bad_value), "object\\$\\.model\\$d_resid_mat cannot contain non-finite values")
})

test_that("print.dnn_diagnostics runs", {
  n <- 30L
  x_mat <- matrix(rnorm(n * 2L), n, 2L)
  y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)
  fit <- np_dnn(
    X = x_mat, x0 = c(0, 0), Y = y_vec,
    method = "exact", k = 2L, s = 3L, se = FALSE, return = "W"
  )
  d <- diagnostics(fit)
  txt <- capture.output(print(d))
  expect_true(any(grepl("<dnn_diagnostics>", txt, fixed = TRUE)))
  expect_true(any(grepl("outcome:", txt, fixed = TRUE)))

  demo_cate <- cate_demo_data(n = 45L, p = 2L, seed = 26L)
  fit_cate <- cate_dnn(
    data = demo_cate$data,
    x_cols = c("X1", "X2"),
    y_col = "Y",
    t_col = "D",
    m_col = "m_oracle",
    e_col = "e_oracle",
    x0 = demo_cate$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )
  txt_cate <- capture.output(print(diagnostics(fit_cate)))
  expect_true(any(grepl("overlap:", txt_cate, fixed = TRUE)))

  demo_plqr <- plqr_demo_data(n = 55L, p = 2L, p_d = 1L, seed = 27L)
  fit_plqr <- plqr_dnn(
    data = demo_plqr$data,
    x_cols = c("X1", "X2"),
    D_cols = "D1",
    q_col = "q_oracle",
    mD_cols = "mD1",
    y_col = "Y",
    beta = demo_plqr$truth$beta_true,
    x0 = demo_plqr$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE,
    return = "W"
  )
  txt_plqr <- capture.output(print(diagnostics(fit_plqr)))
  expect_true(any(grepl("gram   :", txt_plqr, fixed = TRUE)))
})
