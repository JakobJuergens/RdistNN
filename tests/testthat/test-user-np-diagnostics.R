test_that("diagnostics(rdistnn_np_result) returns expected schema and identities", {
  .test_rng_scope({
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

    expect_s3_class(d, "rdistnn_diagnostics")
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
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_np_result) reports unavailable distance fields when dx is absent", {
  .test_rng_scope({
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
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_np_result) errors clearly when W is not retained", {
  .test_rng_scope({
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
      class = "rdistnn_extra_unavailable"
    )
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_np_result) requires current fitted outcome state", {
  .test_rng_scope({
    n <- 35L
    x_mat <- matrix(rnorm(n * 2L), n, 2L)
    y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)

    fit <- np_dnn(
      X = x_mat, x0 = c(0, 0), Y = y_vec,
      method = "exact", k = 2L, s = 3L, se = FALSE,
      return = "W"
    )
    fit$.model <- NULL
    err <- tryCatch(diagnostics(fit), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    expect_identical(err$code, "state_invalid")
    expect_match(err$argument, "object$.model", fixed = TRUE)
    expect_match(conditionMessage(err), "Refit", fixed = TRUE)
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_np_result) validates retained outcome state", {
  .test_rng_scope({
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
    bad_na <- fit
    bad_na$.model$y_vec[1L] <- NA_real_
    for (case in list(list(fit = bad_length, reason = "length"), list(fit = bad_na, reason = "NA"))) {
      error <- tryCatch(diagnostics(case$fit), error = identity)
      expect_s3_class(error, "rdistnn_state_error")
      expect_identical(error$code, "state_invalid")
      expect_identical(error$argument, "object$.model$y_vec")
      expect_match(conditionMessage(error), case$reason, fixed = TRUE)
      expect_match(conditionMessage(error), "[Rr]efit")
    }
  }, seed = 777L)
})
