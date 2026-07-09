skip_if_not(
  all(vapply(
    c("np_dnn", "new_tuning", "tuning_dnn_paper", "rank_weight_dnn_exact"),
    exists, logical(1),
    mode = "function"
  )),
  "np_dnn() grid-mode dependencies not found; skipping file"
)

.layer5_np_data <- function(n = 42L, p = 2L, m = 3L, seed = 5505L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X0 <- matrix(rnorm(m * p), nrow = m, ncol = p)
  y <- 0.75 * X[, 1L] - 0.35 * X[, min(2L, p)] + rnorm(n, sd = 0.2)
  list(X = X, X0 = X0, y = y)
}

.layer5_np_tuning <- function(default_s = 4L, s_grid = c(2L, 4L, 4L, 6L)) {
  new_tuning(
    name = "layer5_np_grid",
    eval = function(...) list(s = default_s, meta = list(lambda = 1)),
    eval_grid = function(...) {
      data.frame(
        lambda = c(0.5, 1, 1.25, 1.5),
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
}

.layer5_top_share <- function(n, k, s, top_n) {
  tab <- rank_weight_dnn_exact(k = k, s = s)$table(n)
  abs_tab <- abs(tab)
  sum(abs_tab[seq_len(min(top_n, length(abs_tab)))]) / sum(abs_tab)
}

.layer5_expect_fixed_top5_diagnostics <- function(grid) {
  diag <- grid$tuning_grid$diagnostics
  expect_true(all(c("top_n", "top5_weight_share", "top_share") %in% names(diag)))
  expect_equal(diag$top_n, rep(5L, nrow(diag)))

  ok_idx <- which(diag$status == "ok")
  expect_gt(length(ok_idx), 0L)
  expect_equal(diag$top5_weight_share[ok_idx], diag$top_share[ok_idx], tolerance = 1e-12)

  differs_from_estimator_k <- logical(length(ok_idx))
  for (pos in seq_along(ok_idx)) {
    i <- ok_idx[[pos]]
    cand <- grid$tuning_grid$candidates[
      grid$tuning_grid$candidates$candidate_id == diag$candidate_id[[i]],
      ,
      drop = FALSE
    ]
    expected_top5 <- .layer5_top_share(
      n = grid$tuning_grid$meta$n,
      k = grid$tuning_grid$meta$k,
      s = cand$s[[1L]],
      top_n = 5L
    )
    estimator_k_share <- .layer5_top_share(
      n = grid$tuning_grid$meta$n,
      k = grid$tuning_grid$meta$k,
      s = cand$s[[1L]],
      top_n = grid$tuning_grid$meta$k
    )
    expect_equal(diag$top5_weight_share[[i]], expected_top5, tolerance = 1e-12)
    differs_from_estimator_k[[pos]] <- !isTRUE(all.equal(expected_top5, estimator_k_share))
  }
  expect_true(any(differs_from_estimator_k))
}

test_that("np_dnn grid mode preserves default-candidate fields and extras", {
  d <- .layer5_np_data()
  tun <- .layer5_np_tuning()

  scalar <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = tun,
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ess", "tuning")
  )
  grid <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ess", "tuning")
  )

  expect_s3_class(grid, "np_result")
  expect_equal(grid$mu, scalar$mu, tolerance = 0)
  expect_equal(grid$se, scalar$se, tolerance = 0)
  expect_equal(grid$W, scalar$W, tolerance = 0)
  expect_equal(grid$ranks, scalar$ranks, tolerance = 0)
  expect_equal(grid$ess, scalar$ess, tolerance = 0)
  expect_equal(grid$meta$s, scalar$meta$s)
  expect_equal(grid$meta$tuning, scalar$meta$tuning)

  expect_named(grid$tuning_grid, c("candidates", "estimates", "diagnostics", "meta"))
  default_id <- grid$tuning_grid$candidates$candidate_id[grid$tuning_grid$candidates$is_default]
  default_rows <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$candidate_id == default_id, ]
  expect_equal(default_rows$estimate, grid$mu, tolerance = 0)
  expect_equal(default_rows$se, grid$se, tolerance = 0)
})

test_that("np_dnn grid mode candidate rows match separate scalar fits", {
  d <- .layer5_np_data(n = 38L, m = 4L, seed = 5506L)
  tun <- .layer5_np_tuning()

  grid <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = character(0)
  )

  ok_candidates <- grid$tuning_grid$candidates[grid$tuning_grid$candidates$status == "ok", ]
  for (i in seq_len(nrow(ok_candidates))) {
    row <- ok_candidates[i, ]
    ref <- np_dnn(
      X = d$X, X0 = d$X0, Y = d$y,
      method = "exact", k = 1L, s = row$s,
      se = TRUE, jackknife = "fast",
      return = character(0)
    )
    got <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$candidate_id == row$candidate_id, ]
    expect_equal(got$estimate, ref$mu, tolerance = 1e-12)
    expect_equal(got$se, ref$se, tolerance = 1e-12)
  }

  dup <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$computed_key == "s_4", ]
  expect_equal(
    dup$estimate[dup$candidate_id == "lambda_1"],
    dup$estimate[dup$candidate_id == "lambda_1p25"],
    tolerance = 0
  )
})

test_that("np_dnn grid diagnostics expose fixed top-5 shares independent of estimator k", {
  d <- .layer5_np_data(n = 40L, m = 3L, seed = 5509L)
  grid <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 2L,
    tuning = .layer5_np_tuning(default_s = 5L, s_grid = c(2L, 5L, 5L, 7L)),
    tuning_mode = "grid",
    se = FALSE,
    return = character(0)
  )

  .layer5_expect_fixed_top5_diagnostics(grid)
})

test_that("np_dnn grid mode computes s = n candidates when no jackknife is requested", {
  n <- 34L
  k <- 2L
  d <- .layer5_np_data(n = n, m = 3L, seed = 5510L)
  grid <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = k,
    tuning = .layer5_np_tuning(default_s = 5L, s_grid = c(2L, 5L, n, n)),
    tuning_mode = "grid",
    se = FALSE,
    return = character(0)
  )

  s_n_ids <- grid$tuning_grid$candidates$candidate_id[grid$tuning_grid$candidates$s == n]
  expect_length(s_n_ids, 2L)
  expect_true(all(grid$tuning_grid$candidates$status[grid$tuning_grid$candidates$s == n] == "ok"))

  s_n_est <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$candidate_id %in% s_n_ids, ]
  expect_true(all(is.finite(s_n_est$estimate)))
  expect_true(all(is.na(s_n_est$se)))

  s_n_diag <- grid$tuning_grid$diagnostics[grid$tuning_grid$diagnostics$candidate_id %in% s_n_ids, ]
  expect_equal(s_n_diag$ess, rep(k, nrow(s_n_diag)), tolerance = 1e-12)
  expect_equal(s_n_diag$max_weight, rep(1 / k, nrow(s_n_diag)), tolerance = 1e-12)
  expect_equal(s_n_diag$top5_weight_share, rep(1, nrow(s_n_diag)), tolerance = 1e-12)
  expect_equal(unique(s_n_diag$computed_key), sprintf("s_%d", n))
})

test_that("np_dnn grid mode rejects explicit s and can use the paper default grid", {
  d <- .layer5_np_data(n = 40L, m = 2L, seed = 5507L)

  expect_error(
    np_dnn(
      X = d$X, X0 = d$X0, Y = d$y,
      method = "exact", k = 1L, s = 4L,
      tuning_mode = "grid", se = FALSE
    ),
    "tuning_mode.*grid.*s"
  )

  grid <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 1L,
    tuning_mode = "grid", se = FALSE,
    return = character(0)
  )

  expect_true(is.list(grid$tuning_grid))
  expect_equal(grid$tuning_grid$candidates$lambda, c(0.75, 1, 1.25, 1.5))
  expect_equal(grid$tuning_grid$meta$name, "paper_undersmooth")
  expect_equal(grid$meta$tuning$source, "paper_default")
})

test_that("np_dnn grid mode records infeasible candidates with one warning", {
  d <- .layer5_np_data(n = 36L, m = 2L, seed = 5508L)
  tun <- .layer5_np_tuning(default_s = 4L, s_grid = c(1L, 4L, 4L, 36L))

  expect_warning(
    grid <- np_dnn(
      X = d$X, X0 = d$X0, Y = d$y,
      method = "exact", k = 2L, tuning = tun, tuning_mode = "grid",
      se = TRUE, jackknife = "fast",
      return = character(0)
    ),
    "Skipped infeasible tuning candidate"
  )

  candidates <- grid$tuning_grid$candidates
  expect_equal(candidates$status, c("infeasible", "ok", "ok", "infeasible"))
  skipped <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$status == "infeasible", ]
  expect_true(all(is.na(skipped$estimate)))
  expect_true(all(is.na(skipped$se)))
  expect_equal(
    grid$tuning_grid$diagnostics$status,
    c("infeasible", "ok", "ok", "infeasible")
  )
})

test_that("later layers add tuning_mode to all main DNN wrappers", {
  expect_true("tuning_mode" %in% names(formals(np_dnn)))
  expect_true("tuning_mode" %in% names(formals(cate_dnn)))
  expect_true("tuning_mode" %in% names(formals(plqr_dnn)))
})
