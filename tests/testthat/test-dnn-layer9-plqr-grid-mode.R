skip_if_not(
  all(vapply(
    c(
      "plqr_dnn", "new_tuning", "tuning_dnn_paper",
      "candidate_fit", "candidate_fits", "rank_weight_dnn_exact"
    ),
    exists, logical(1),
    mode = "function"
  )),
  "plqr_dnn() grid-mode dependencies not found; skipping file"
)

.layer9_plqr_tuning <- function(default_s = 5L, s_grid = c(3L, 5L, 5L, 7L)) {
  new_tuning(
    name = "layer9_plqr_grid",
    eval = function(...) list(s = default_s, meta = list(kappa = 1)),
    eval_grid = function(...) {
      data.frame(
        kappa = c(0.6, 1, 1.25, 1.5),
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
}

.layer9_plqr_case <- function(n = 46L, p = 2L, p_d = 2L, m = 3L, seed = 9909L) {
  demo <- plqr_demo_data(n = n, p = p, p_d = p_d, tau = 0.5, seed = seed)
  df <- demo$data
  X0 <- matrix(
    seq(-0.6, 0.6, length.out = m * p),
    nrow = m,
    ncol = p
  )
  list(
    df = df,
    X = as.matrix(df[paste0("X", seq_len(p))]),
    D = as.matrix(df[paste0("D", seq_len(p_d))]),
    q = as.numeric(df$q_oracle),
    mD = as.matrix(df[paste0("mD", seq_len(p_d))]),
    Y = as.numeric(df$Y),
    X0 = X0,
    beta = demo$truth$beta_true,
    tau = 0.5,
    x_cols = paste0("X", seq_len(p)),
    D_cols = paste0("D", seq_len(p_d)),
    mD_cols = paste0("mD", seq_len(p_d))
  )
}

.layer9_do_plqr <- function(d, ...) {
  plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    ...
  )
}

.layer9_top_share <- function(n, k, s, top_n) {
  tab <- rank_weight_dnn_exact(k = k, s = s)$table(n)
  abs_tab <- abs(tab)
  sum(abs_tab[seq_len(min(top_n, length(abs_tab)))]) / sum(abs_tab)
}

.layer9_expect_fixed_top5_diagnostics <- function(grid) {
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
    expected_top5 <- .layer9_top_share(
      n = grid$tuning_grid$meta$n,
      k = grid$tuning_grid$meta$k,
      s = cand$s[[1L]],
      top_n = 5L
    )
    estimator_k_share <- .layer9_top_share(
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

.layer9_value_matrix <- function(rows, value_col, p_d, m) {
  comps <- paste0("comp_", seq_len(p_d))
  out <- matrix(NA_real_, nrow = p_d, ncol = m)
  rownames(out) <- comps
  colnames(out) <- paste0("target_", seq_len(m))
  for (i in seq_len(nrow(rows))) {
    out[match(rows$component[[i]], comps), as.integer(rows$target_id[[i]])] <-
      rows[[value_col]][[i]]
  }
  out
}

.layer9_expect_plqr_candidate_rows_match <- function(rows, ref, candidate_id) {
  got_moment <- .layer9_value_matrix(
    rows,
    "estimate",
    p_d = nrow(ref$moment),
    m = ncol(ref$moment)
  )
  got_se <- .layer9_value_matrix(
    rows,
    "se",
    p_d = nrow(ref$moment),
    m = ncol(ref$moment)
  )
  expect_equal(got_moment, ref$moment, tolerance = 1e-12, info = candidate_id)
  expect_equal(got_se, ref$se, tolerance = 1e-12, info = candidate_id)
}

test_that("plqr_dnn grid mode preserves default-candidate fields and long-form rows", {
  d <- .layer9_plqr_case(seed = 9910L)
  tun <- .layer9_plqr_tuning()

  scalar <- .layer9_do_plqr(
    d,
    method = "exact", k = 1L, tuning = tun,
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ess", "g", "tuning")
  )
  grid <- .layer9_do_plqr(
    d,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ess", "g", "tuning")
  )

  expect_s3_class(grid, "plqr_result")
  expect_equal(grid$moment, scalar$moment, tolerance = 0)
  expect_equal(grid$se, scalar$se, tolerance = 0)
  expect_equal(grid$W, scalar$W, tolerance = 0)
  expect_equal(grid$ranks, scalar$ranks, tolerance = 0)
  expect_equal(grid$ess, scalar$ess, tolerance = 0)
  expect_equal(grid$g, scalar$g, tolerance = 0)
  expect_equal(grid$meta$s, scalar$meta$s)
  expect_equal(grid$meta$tuning, scalar$meta$tuning)

  expect_named(grid$tuning_grid, c("candidates", "estimates", "diagnostics", "meta"))
  expect_equal(
    nrow(grid$tuning_grid$estimates),
    nrow(grid$tuning_grid$candidates) * nrow(grid$moment) * ncol(grid$moment)
  )
  expect_equal(sort(unique(grid$tuning_grid$estimates$component)), rownames(grid$moment))

  default_id <- grid$tuning_grid$candidates$candidate_id[grid$tuning_grid$candidates$is_default]
  default_rows <- grid$tuning_grid$estimates[
    grid$tuning_grid$estimates$candidate_id == default_id,
  ]
  .layer9_expect_plqr_candidate_rows_match(default_rows, grid, default_id)
})

test_that("plqr_dnn grid mode candidate rows match separate scalar fits", {
  d <- .layer9_plqr_case(n = 44L, m = 4L, seed = 9911L)
  tun <- .layer9_plqr_tuning()
  grid <- .layer9_do_plqr(
    d,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = character(0)
  )

  ok_candidates <- grid$tuning_grid$candidates[grid$tuning_grid$candidates$status == "ok", ]
  for (i in seq_len(nrow(ok_candidates))) {
    row <- ok_candidates[i, ]
    ref <- .layer9_do_plqr(
      d,
      method = "exact", k = 1L, s = row$s,
      se = TRUE, jackknife = "fast",
      return = character(0)
    )
    got <- grid$tuning_grid$estimates[
      grid$tuning_grid$estimates$candidate_id == row$candidate_id,
    ]
    .layer9_expect_plqr_candidate_rows_match(got, ref, row$candidate_id)
  }

  dup <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$computed_key == "s_5", ]
  expect_equal(
    dup$estimate[dup$candidate_id == "kappa_1"],
    dup$estimate[dup$candidate_id == "kappa_1p25"],
    tolerance = 0
  )
  expect_equal(
    dup$se[dup$candidate_id == "kappa_1"],
    dup$se[dup$candidate_id == "kappa_1p25"],
    tolerance = 0
  )
})

test_that("plqr_dnn grid diagnostics expose fixed top-5 shares independent of estimator k", {
  d <- .layer9_plqr_case(n = 44L, m = 3L, seed = 9914L)
  grid <- .layer9_do_plqr(
    d,
    method = "exact", k = 2L,
    tuning = .layer9_plqr_tuning(default_s = 5L, s_grid = c(2L, 5L, 5L, 7L)),
    tuning_mode = "grid",
    se = FALSE,
    return = character(0)
  )

  .layer9_expect_fixed_top5_diagnostics(grid)
})

test_that("candidate_fit reconstructs PLQR candidate views", {
  d <- .layer9_plqr_case(n = 42L, seed = 9912L)
  tun <- .layer9_plqr_tuning()
  fit <- .layer9_do_plqr(
    d,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "g")
  )

  cand <- fit$tuning_grid$candidates[fit$tuning_grid$candidates$candidate_id == "kappa_1p5", ]
  view <- candidate_fit(fit, cand$candidate_id)
  ref <- .layer9_do_plqr(
    d,
    method = "exact", k = 1L, s = cand$s,
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "ess", "g")
  )

  expect_s3_class(view, "plqr_result")
  expect_equal(view$moment, ref$moment, tolerance = 1e-12)
  expect_equal(view$se, ref$se, tolerance = 1e-12)
  expect_equal(view$W, ref$W, tolerance = 1e-12)
  expect_equal(view$ess, ref$ess, tolerance = 1e-12)
  expect_equal(view$ranks, ref$ranks, tolerance = 0)
  expect_equal(view$ord, ref$ord, tolerance = 0)
  expect_equal(view$dx, ref$dx, tolerance = 1e-12)
  expect_equal(view$g, ref$g, tolerance = 0)
  expect_equal(view$meta$s, ref$meta$s)
  expect_equal(view$meta$candidate_id, cand$candidate_id)
  expect_null(view$tuning_grid)

  views <- candidate_fits(fit, candidates = c("kappa_1p5", "kappa_1"))
  expect_equal(names(views), c("kappa_1p5", "kappa_1"))
  expect_equal(views[[1L]], view)
})

test_that("plqr_dnn grid mode rejects explicit s and records infeasible candidates", {
  d <- .layer9_plqr_case(n = 36L, m = 2L, seed = 9913L)

  expect_error(
    .layer9_do_plqr(
      d,
      method = "exact", k = 1L, s = 4L,
      tuning_mode = "grid", se = FALSE
    ),
    "tuning_mode.*grid.*s"
  )

  grid <- .layer9_do_plqr(
    d,
    method = "exact", k = 1L,
    tuning_mode = "grid", se = FALSE,
    return = character(0)
  )
  expect_true(is.list(grid$tuning_grid))
  expect_equal(grid$tuning_grid$candidates$kappa, c(0.75, 1, 1.25, 1.5))
  expect_equal(grid$tuning_grid$meta$name, "paper_undersmooth")
  expect_equal(grid$meta$tuning$source, "paper_default")

  bad_tun <- .layer9_plqr_tuning(default_s = 5L, s_grid = c(1L, 5L, 5L, 36L))
  expect_warning(
    bad <- .layer9_do_plqr(
      d,
      method = "exact", k = 2L, tuning = bad_tun, tuning_mode = "grid",
      se = TRUE, jackknife = "fast",
      return = character(0)
    ),
    "Skipped infeasible tuning candidate"
  )
  expect_equal(bad$tuning_grid$candidates$status, c("infeasible", "ok", "ok", "infeasible"))
  skipped <- bad$tuning_grid$estimates[bad$tuning_grid$estimates$status == "infeasible", ]
  expect_true(all(is.na(skipped$estimate)))
  expect_true(all(is.na(skipped$se)))
})

test_that("Layer 9 adds tuning_mode to all main DNN wrappers", {
  expect_true("tuning_mode" %in% names(formals(np_dnn)))
  expect_true("tuning_mode" %in% names(formals(cate_dnn)))
  expect_true("tuning_mode" %in% names(formals(plqr_dnn)))
})
