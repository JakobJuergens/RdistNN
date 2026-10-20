.candidate_view_tuning <- function(default_s = 4L, s_grid = c(2L, 4L, 6L)) {
  new_tuning(
    name = "layer7_grid",
    eval = function(...) list(s = default_s, meta = list(lambda = 1)),
    eval_grid = function(...) {
      data.frame(
        lambda = c(0.5, 1, 1.5),
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
}

.candidate_view_np_case <- function(n = 42L, p = 2L, m = 3L, seed = 7707L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X0 <- matrix(rnorm(m * p), nrow = m, ncol = p)
  y <- X[, 1L] - 0.35 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
  list(X = X, X0 = X0, y = y)
}

.candidate_view_rank_matrix <- function(n, m) {
  matrix(rep(seq_len(n), m), nrow = n, ncol = m)
}

.candidate_view_compare_np <- function(view, ref) {
  expect_s3_class(view, "np_result")
  expect_equal(view$mu, ref$mu, tolerance = 1e-12)
  expect_equal(view$se, ref$se, tolerance = 1e-12)
  expect_equal(view$W, ref$W, tolerance = 1e-12)
  expect_equal(view$ess, ref$ess, tolerance = 1e-12)
  expect_equal(view$ranks, ref$ranks, tolerance = 0)
  expect_equal(view$ord, ref$ord, tolerance = 0)
  expect_equal(view$dx, ref$dx, tolerance = 1e-12)
  expect_equal(view$y, ref$y, tolerance = 0)
  expect_equal(view$meta$s, ref$meta$s)
}

.candidate_view_compare_cate <- function(view, ref) {
  expect_s3_class(view, "cate_result")
  expect_equal(view$tau, ref$tau, tolerance = 1e-12)
  expect_equal(view$se, ref$se, tolerance = 1e-12)
  expect_equal(view$W, ref$W, tolerance = 1e-12)
  expect_equal(view$ess, ref$ess, tolerance = 1e-12)
  expect_equal(view$ranks, ref$ranks, tolerance = 0)
  expect_equal(view$ord, ref$ord, tolerance = 0)
  expect_equal(view$dx, ref$dx, tolerance = 1e-12)
  expect_equal(view$psi, ref$psi, tolerance = 0)
  expect_equal(view$nuisance, ref$nuisance, tolerance = 0)
  expect_equal(view$meta$s, ref$meta$s)
}

test_that("candidate_fit reconstructs NP candidate views equivalent to scalar fits", {
  d <- .candidate_view_np_case(seed = 7710L)
  tun <- .candidate_view_tuning()
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "neighbors", "y")
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "lambda_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, s = cand$s,
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "neighbors", "y")
  )

  .candidate_view_compare_np(view, ref)
  expect_equal(view$meta$candidate_id, cand$candidate_id)
  expect_equal(view$neighbors, ref$neighbors)
  expect_null(view$tuning_grid)

  default_id <- fit$tuning_grid$candidates$candidate_id[fit$tuning_grid$candidates$is_default]
  default_view <- candidate_fit(fit, default_id)
  expect_equal(default_view$mu, fit$mu, tolerance = 0)
  expect_equal(default_view$se, fit$se, tolerance = 0)
  expect_equal(default_view$W, fit$W, tolerance = 0)
})

test_that("candidate_fit preserves matrix extras for single-target x0 views", {
  d <- .candidate_view_np_case(n = 35L, m = 1L, seed = 7715L)
  x0 <- as.numeric(d$X0[1L, ])
  tun <- .candidate_view_tuning(default_s = 4L, s_grid = c(2L, 4L, 7L))
  fit <- np_dnn(
    X = d$X, x0 = x0, Y = d$y,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = FALSE,
    return = c("W", "ranks", "ord", "dx", "ess", "neighbors", "y")
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "lambda_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- np_dnn(
    X = d$X, x0 = x0, Y = d$y,
    method = "exact", k = 1L, s = cand$s,
    se = FALSE,
    return = c("W", "ranks", "ord", "dx", "ess", "neighbors", "y")
  )

  .candidate_view_compare_np(view, ref)
  expect_equal(dim(view$x0), c(1L, ncol(d$X)))
  expect_equal(dim(view$W), c(nrow(d$X), 1L))
  expect_equal(dim(view$ranks), c(nrow(d$X), 1L))
  expect_equal(dim(view$ord), c(nrow(d$X), 1L))
  expect_equal(dim(view$dx), c(nrow(d$X), 1L))
  expect_equal(dim(view$neighbors), c(1L, 1L))
  expect_length(view$ess, 1L)
  expect_length(view$mu, 1L)
  expect_equal(view$neighbors[1L, 1L], order(view$ranks[, 1L])[[1L]])
  expect_equal(sum(view$W[, 1L] > 0), nrow(d$X) - cand$s[[1L]] + 1L)
})

test_that("candidate_fit keeps nearest-k neighbors distinct from DNN support", {
  d <- .candidate_view_np_case(n = 42L, m = 4L, seed = 7714L)
  k <- 2L
  tun <- .candidate_view_tuning(default_s = 5L, s_grid = c(3L, 5L, 8L))
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = k, tuning = tun, tuning_mode = "grid",
    se = FALSE,
    return = c("W", "ranks", "neighbors")
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "lambda_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = k, s = cand$s,
    se = FALSE,
    return = c("W", "ranks", "neighbors")
  )

  support_edge <- nrow(d$X) - cand$s[[1L]] + k
  expect_equal(view$W, ref$W, tolerance = 1e-12)
  expect_equal(view$neighbors, ref$neighbors, tolerance = 0)
  expect_equal(dim(view$neighbors), c(k, nrow(d$X0)))
  for (j in seq_len(nrow(d$X0))) {
    expect_equal(view$neighbors[, j], order(view$ranks[, j])[seq_len(k)])
    expect_equal(sum(view$W[, j] > 0), support_edge)
    expect_gt(sum(view$W[, j] > 0), nrow(view$neighbors))
  }
})

test_that("candidate_fit reconstructs CATE candidate views and methods work", {
  d <- .make_df_data(n = 44L, p = 2L, m = 3L, seed = 7711L)
  tun <- .candidate_view_tuning()
  fit <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
    X0 = d$X0,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "psi", "nuisance")
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "lambda_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
    X0 = d$X0,
    method = "exact", k = 1L, s = cand$s,
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "psi", "nuisance")
  )

  .candidate_view_compare_cate(view, ref)
  expect_equal(view$meta$candidate_id, cand$candidate_id)
  expect_null(view$tuning_grid)

  expect_equal(predict(view, X0 = d$X0), predict(ref, X0 = d$X0), tolerance = 1e-12)
  diag_view <- diagnostics(view)
  diag_ref <- diagnostics(ref)
  expect_equal(diag_view$n_eff, diag_ref$n_eff, tolerance = 1e-12)
  expect_equal(diag_view$s, rep(cand$s, nrow(d$X0)))
  expect_equal(diag_view$overlap_available, diag_ref$overlap_available)
  expect_equal(diag_view$propensity_available, diag_ref$propensity_available)
})

test_that("CATE candidate prediction works without retained W but diagnostics require it", {
  d <- .make_df_data(n = 44L, p = 2L, m = 3L, seed = 7717L)
  tun <- .candidate_view_tuning()
  fit <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
    X0 = d$X0,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = FALSE,
    return = character(0)
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "lambda_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
    X0 = d$X0,
    method = "exact", k = 1L, s = cand$s,
    se = FALSE,
    return = character(0)
  )

  expect_false("W" %in% names(view))
  expect_false("ess" %in% names(view))
  expect_equal(predict(view, X0 = d$X0), view$tau, tolerance = 1e-12)
  expect_equal(predict(view, X0 = d$X0), predict(ref, X0 = d$X0), tolerance = 1e-12)
  expect_error(
    diagnostics(view),
    "refit with `return` including \"W\""
  )
})

test_that("candidate_fits is a named list matching repeated candidate_fit calls", {
  d <- .candidate_view_np_case(seed = 7712L)
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = .candidate_view_tuning(), tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ess")
  )

  ids <- fit$tuning_grid$candidates$candidate_id
  views <- candidate_fits(fit)
  expect_type(views, "list")
  expect_equal(names(views), ids)
  expect_equal(views[[ids[[1L]]]], candidate_fit(fit, ids[[1L]]))

  subset_views <- candidate_fits(fit, candidates = ids[c(3L, 1L)])
  expect_equal(names(subset_views), ids[c(3L, 1L)])
  expect_equal(subset_views[[1L]], candidate_fit(fit, ids[[3L]]))
})

test_that("candidate_fits rejects empty, malformed, and unknown selections", {
  d <- .candidate_view_np_case(seed = 7718L)
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = .candidate_view_tuning(), tuning_mode = "grid",
    se = FALSE,
    return = character(0)
  )
  ids <- fit$tuning_grid$candidates$candidate_id

  expect_error(
    candidate_fits(fit, candidates = character(0)),
    "one or more non-missing candidate id"
  )
  expect_error(
    candidate_fits(fit, candidates = c(ids[[1L]], "")),
    "one or more non-missing candidate id"
  )
  expect_error(
    candidate_fits(fit, candidates = c(ids[[1L]], NA_character_)),
    "one or more non-missing candidate id"
  )
  expect_error(
    candidate_fits(fit, candidates = "not_a_candidate"),
    "unknown candidate id"
  )
  expect_error(
    candidate_fits(fit, candidates = numeric(0)),
    "one or more integer row indices"
  )
})

test_that("candidate_fit does not attach optional extras absent from the parent fit", {
  d <- .candidate_view_np_case(seed = 7716L)
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L,
    tuning = .candidate_view_tuning(), tuning_mode = "grid",
    se = FALSE,
    return = character(0)
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "lambda_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, s = cand$s,
    se = FALSE,
    return = character(0)
  )

  expect_equal(view$mu, ref$mu, tolerance = 1e-12)
  expect_null(view$se)
  expect_equal(view$meta$s, cand$s[[1L]])
  expect_equal(view$meta$candidate_id, cand$candidate_id)
  expect_null(view$tuning_grid)
  expect_false(any(c("W", "ranks", "ord", "dx", "ess", "neighbors", "y") %in% names(view)))
})

test_that("candidate_fit rejects stale cached ranks inconsistent with result dimensions", {
  d <- .candidate_view_np_case(n = 36L, m = 3L, seed = 7719L)
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L,
    tuning = .candidate_view_tuning(), tuning_mode = "grid",
    se = FALSE,
    return = c("W", "ranks", "ess")
  )

  short_ranks <- fit
  short_ranks$.model$ranks <- .candidate_view_rank_matrix(nrow(d$X) - 1L, nrow(d$X0))
  expect_error(
    candidate_fit(short_ranks, "lambda_1p5"),
    "row count.*object\\$meta\\$n"
  )

  few_target_ranks <- fit
  few_target_ranks$.model$ranks <- .candidate_view_rank_matrix(nrow(d$X), nrow(d$X0) - 1L)
  expect_error(
    candidate_fit(few_target_ranks, "lambda_1p5"),
    "column count.*target points"
  )

  stale_x0 <- fit
  stale_x0$x0 <- stale_x0$x0[-1L, , drop = FALSE]
  expect_error(
    candidate_fit(stale_x0, "lambda_1p5"),
    "object\\$x0.*object\\$meta\\$m"
  )

  fallback_ranks <- fit
  fallback_ranks$.model$ranks <- NULL
  fallback_ranks$ranks <- .candidate_view_rank_matrix(nrow(d$X) - 1L, nrow(d$X0))
  expect_error(
    candidate_fit(fallback_ranks, "lambda_1p5"),
    "row count.*object\\$meta\\$n"
  )
})

test_that("candidate_fit rejects stale candidate estimate rows", {
  d <- .candidate_view_np_case(n = 36L, m = 3L, seed = 7720L)
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L,
    tuning = .candidate_view_tuning(), tuning_mode = "grid",
    se = FALSE,
    return = character(0)
  )

  missing_target <- fit
  drop <- missing_target$tuning_grid$estimates$candidate_id == "lambda_1p5" &
    missing_target$tuning_grid$estimates$target_id == 3L
  missing_target$tuning_grid$estimates <- missing_target$tuning_grid$estimates[!drop, , drop = FALSE]
  expect_error(
    candidate_fit(missing_target, "lambda_1p5"),
    "missing target_id"
  )

  duplicate_target <- fit
  dup <- duplicate_target$tuning_grid$estimates[
    duplicate_target$tuning_grid$estimates$candidate_id == "lambda_1p5" &
      duplicate_target$tuning_grid$estimates$target_id == 1L,
    ,
    drop = FALSE
  ]
  duplicate_target$tuning_grid$estimates <- rbind(duplicate_target$tuning_grid$estimates, dup)
  expect_error(
    candidate_fit(duplicate_target, "lambda_1p5"),
    "duplicate target_id"
  )

  invalid_target <- fit
  invalid_idx <- which(
    invalid_target$tuning_grid$estimates$candidate_id == "lambda_1p5"
  )[[1L]]
  invalid_target$tuning_grid$estimates$target_id[[invalid_idx]] <- 99L
  expect_error(
    candidate_fit(invalid_target, "lambda_1p5"),
    "target_id.*outside 1..3"
  )
})

test_that("candidate_fit errors are actionable", {
  d <- .candidate_view_np_case(n = 30L, m = 2L, seed = 7713L)
  scalar <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, s = 4L,
    se = FALSE
  )
  expect_error(
    candidate_fit(scalar, "lambda_1"),
    "tuning_mode = \"grid\""
  )

  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = .candidate_view_tuning(), tuning_mode = "grid",
    se = FALSE
  )
  expect_error(candidate_fit(fit, "not_a_candidate"), "unknown candidate")
  expect_error(candidate_fit(fit, c("lambda_1", "lambda_1p5")), "single")
  expect_error(candidate_fit(fit, 0L), "positive integer|row index")
  expect_error(candidate_fit(fit, 99L), "row index must be in 1..")

  missing_est <- fit
  missing_est$tuning_grid$estimates <- missing_est$tuning_grid$estimates[
    missing_est$tuning_grid$estimates$candidate_id != "lambda_1p5",
    ,
    drop = FALSE
  ]
  expect_error(candidate_fit(missing_est, "lambda_1p5"), "No estimates recorded")

  malformed <- fit
  malformed$tuning_grid$diagnostics <- NULL
  expect_error(candidate_fit(malformed, "lambda_1"), "tuning_grid.*malformed")

  bad_tun <- .candidate_view_tuning(default_s = 4L, s_grid = c(1L, 4L, 6L))
  infeasible <- suppressWarnings(np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 2L, tuning = bad_tun, tuning_mode = "grid",
    se = TRUE,
    return = "W"
  ))
  bad_id <- infeasible$tuning_grid$candidates$candidate_id[
    infeasible$tuning_grid$candidates$status == "infeasible"
  ][[1L]]
  expect_error(candidate_fit(infeasible, bad_id), "infeasible.*k > s")

  failed <- fit
  failed$tuning_grid$candidates$status[[1L]] <- "failed"
  failed$tuning_grid$candidates$reason[[1L]] <- "boom in candidate"
  expect_error(candidate_fit(failed, 1L), "failed.*boom in candidate")
})
