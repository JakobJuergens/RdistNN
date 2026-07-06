skip_if_not(
  all(vapply(
    c("np_dnn", "cate_dnn", "new_tuning", "tuning_dnn_paper", "diagnostics"),
    exists, logical(1),
    mode = "function"
  )),
  "Layer 7 candidate-fit dependencies not found; skipping file"
)

.layer7_tuning <- function(default_s = 4L, s_grid = c(2L, 4L, 6L)) {
  new_tuning(
    name = "layer7_grid",
    eval = function(...) list(s = default_s, meta = list(kappa = 1)),
    eval_grid = function(...) {
      data.frame(
        kappa = c(0.5, 1, 1.5),
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
}

.layer7_np_case <- function(n = 42L, p = 2L, m = 3L, seed = 7707L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X0 <- matrix(rnorm(m * p), nrow = m, ncol = p)
  y <- X[, 1L] - 0.35 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
  list(X = X, X0 = X0, y = y)
}

.layer7_compare_np <- function(view, ref) {
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

.layer7_compare_cate <- function(view, ref) {
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
  d <- .layer7_np_case(seed = 7710L)
  tun <- .layer7_tuning()
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "neighbors", "y")
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "kappa_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, s = cand$s,
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "neighbors", "y")
  )

  .layer7_compare_np(view, ref)
  expect_equal(view$meta$candidate_id, cand$candidate_id)
  expect_equal(view$neighbors, ref$neighbors)
  expect_null(view$tuning_grid)

  default_id <- fit$tuning_grid$candidates$candidate_id[fit$tuning_grid$candidates$is_default]
  default_view <- candidate_fit(fit, default_id)
  expect_equal(default_view$mu, fit$mu, tolerance = 0)
  expect_equal(default_view$se, fit$se, tolerance = 0)
  expect_equal(default_view$W, fit$W, tolerance = 0)
})

test_that("candidate_fit reconstructs CATE candidate views and methods work", {
  d <- .make_df_data(n = 44L, p = 2L, m = 3L, seed = 7711L)
  tun <- .layer7_tuning()
  fit <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m_col = "m_oracle", e_col = "e_oracle",
    X0 = d$X0,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "psi", "nuisance")
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "kappa_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m_col = "m_oracle", e_col = "e_oracle",
    X0 = d$X0,
    method = "exact", k = 1L, s = cand$s,
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "psi", "nuisance")
  )

  .layer7_compare_cate(view, ref)
  expect_equal(view$meta$candidate_id, cand$candidate_id)
  expect_null(view$tuning_grid)

  expect_equal(predict(view, X0 = d$X0), predict(ref, X0 = d$X0), tolerance = 1e-12)
  diag_view <- diagnostics(view)
  diag_ref <- diagnostics(ref)
  expect_equal(diag_view$n_eff, diag_ref$n_eff, tolerance = 1e-12)
  expect_equal(diag_view$s, rep(cand$s, nrow(d$X0)))
  expect_equal(diag_view$treatment_overlap_available, diag_ref$treatment_overlap_available)
})

test_that("candidate_fits is a named list matching repeated candidate_fit calls", {
  d <- .layer7_np_case(seed = 7712L)
  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = .layer7_tuning(), tuning_mode = "grid",
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

test_that("candidate_fit errors are actionable", {
  d <- .layer7_np_case(n = 30L, m = 2L, seed = 7713L)
  scalar <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, s = 4L,
    se = FALSE
  )
  expect_error(
    candidate_fit(scalar, "kappa_1"),
    "tuning_mode = \"grid\""
  )

  fit <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = .layer7_tuning(), tuning_mode = "grid",
    se = FALSE
  )
  expect_error(candidate_fit(fit, "not_a_candidate"), "unknown candidate")
  expect_error(candidate_fit(fit, c("kappa_1", "kappa_1p5")), "single")

  bad_tun <- .layer7_tuning(default_s = 4L, s_grid = c(1L, 4L, 6L))
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
