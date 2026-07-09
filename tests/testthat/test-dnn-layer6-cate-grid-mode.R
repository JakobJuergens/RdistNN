skip_if_not(
  all(vapply(
    c("cate_dnn", "new_tuning", "tuning_dnn_paper", "rank_weight_dnn_exact"),
    exists, logical(1),
    mode = "function"
  )),
  "cate_dnn() grid-mode dependencies not found; skipping file"
)

.layer6_cate_tuning <- function(default_s = 4L, s_grid = c(2L, 4L, 4L, 6L)) {
  new_tuning(
    name = "layer6_cate_grid",
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

.layer6_cate_case <- function(kind = c("direct", "aipw", "robinson"),
                              n = 42L, p = 2L, m = 3L, seed = 6606L) {
  kind <- match.arg(kind)
  d <- .make_df_data(n = n, p = p, m = m, seed = seed)
  x_mat <- as.matrix(d$df[, d$x_cols, drop = FALSE])
  if (identical(kind, "direct")) {
    return(list(
      grid_args = list(X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw),
      scalar_args = list(X = x_mat, X0 = d$X0, pseudo_outcome = d$psi_aipw)
    ))
  }
  if (identical(kind, "aipw")) {
    args <- list(
      data = d$df, x_cols = d$x_cols,
      y_col = "Y", t_col = "D",
      m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
      X0 = d$X0
    )
    return(list(grid_args = args, scalar_args = args))
  }
  args <- list(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D",
    m_col = "m_oracle", e_col = "e_oracle",
    X0 = d$X0
  )
  list(grid_args = args, scalar_args = args)
}

.layer6_do_cate <- function(args, ...) {
  do.call(cate_dnn, c(args, list(...)))
}

.layer6_top_share <- function(n, k, s, top_n) {
  tab <- rank_weight_dnn_exact(k = k, s = s)$table(n)
  abs_tab <- abs(tab)
  sum(abs_tab[seq_len(min(top_n, length(abs_tab)))]) / sum(abs_tab)
}

.layer6_expect_fixed_top5_diagnostics <- function(grid) {
  diag <- grid$tuning_grid$diagnostics
  expect_true(all(c("top_n", "top5_weight_share") %in% names(diag)))
  expect_false("top_share" %in% names(diag))
  expect_equal(diag$top_n, rep(5L, nrow(diag)))

  ok_idx <- which(diag$status == "ok")
  expect_gt(length(ok_idx), 0L)

  differs_from_estimator_k <- logical(length(ok_idx))
  for (pos in seq_along(ok_idx)) {
    i <- ok_idx[[pos]]
    cand <- grid$tuning_grid$candidates[
      grid$tuning_grid$candidates$candidate_id == diag$candidate_id[[i]],
      ,
      drop = FALSE
    ]
    expected_top5 <- .layer6_top_share(
      n = grid$tuning_grid$meta$n,
      k = grid$tuning_grid$meta$k,
      s = cand$s[[1L]],
      top_n = 5L
    )
    estimator_k_share <- .layer6_top_share(
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

.layer6_expect_grid_matches_scalar <- function(kind, seed) {
  case <- .layer6_cate_case(kind = kind, seed = seed)
  tun <- .layer6_cate_tuning()
  grid <- .layer6_do_cate(
    case$grid_args,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = character(0)
  )

  ok_candidates <- grid$tuning_grid$candidates[grid$tuning_grid$candidates$status == "ok", ]
  for (i in seq_len(nrow(ok_candidates))) {
    row <- ok_candidates[i, ]
    ref <- .layer6_do_cate(
      case$scalar_args,
      method = "exact", k = 1L, s = row$s,
      se = TRUE, jackknife = "fast",
      return = character(0)
    )
    got <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$candidate_id == row$candidate_id, ]
    expect_equal(got$estimate, ref$tau, tolerance = 1e-12)
    expect_equal(got$se, ref$se, tolerance = 1e-12)
  }

  dup <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$computed_key == "s_4", ]
  expect_equal(
    dup$estimate[dup$candidate_id == "lambda_1"],
    dup$estimate[dup$candidate_id == "lambda_1p25"],
    tolerance = 0
  )
  invisible(grid)
}

.layer6_count_traced_calls <- function(fun, expr) {
  trace_env <- new.env(parent = emptyenv())
  trace_env$calls <- 0L
  assign(".layer6_trace_env", trace_env, envir = .GlobalEnv)
  on.exit(rm(".layer6_trace_env", envir = .GlobalEnv), add = TRUE)
  suppressMessages(invisible(trace(
    fun,
    tracer = quote({
      .GlobalEnv$.layer6_trace_env$calls <- .GlobalEnv$.layer6_trace_env$calls + 1L
    }),
    where = asNamespace("RdistNN"),
    print = FALSE
  )))
  on.exit(suppressMessages(untrace(fun, where = asNamespace("RdistNN"))), add = TRUE)
  force(expr)
  trace_env$calls
}

test_that("cate_dnn grid mode preserves default-candidate fields and extras", {
  case <- .layer6_cate_case("direct", seed = 6610L)
  tun <- .layer6_cate_tuning()

  scalar <- .layer6_do_cate(
    case$scalar_args,
    method = "exact", k = 1L, tuning = tun,
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ess", "psi", "tuning")
  )
  grid <- .layer6_do_cate(
    case$grid_args,
    method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
    se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "ess", "psi", "tuning")
  )

  expect_s3_class(grid, "cate_result")
  expect_equal(grid$tau, scalar$tau, tolerance = 0)
  expect_equal(grid$se, scalar$se, tolerance = 0)
  expect_equal(grid$W, scalar$W, tolerance = 0)
  expect_equal(grid$ranks, scalar$ranks, tolerance = 0)
  expect_equal(grid$ess, scalar$ess, tolerance = 0)
  expect_equal(grid$psi, scalar$psi, tolerance = 0)
  expect_equal(grid$meta$s, scalar$meta$s)
  expect_equal(grid$meta$tuning, scalar$meta$tuning)

  expect_named(grid$tuning_grid, c("candidates", "estimates", "diagnostics", "meta"))
  default_id <- grid$tuning_grid$candidates$candidate_id[grid$tuning_grid$candidates$is_default]
  default_rows <- grid$tuning_grid$estimates[grid$tuning_grid$estimates$candidate_id == default_id, ]
  expect_equal(default_rows$estimate, grid$tau, tolerance = 0)
  expect_equal(default_rows$se, grid$se, tolerance = 0)
})

test_that("cate_dnn grid mode candidate rows match scalar fits on every pseudo-outcome path", {
  .layer6_expect_grid_matches_scalar("direct", seed = 6611L)
  .layer6_expect_grid_matches_scalar("aipw", seed = 6612L)
  .layer6_expect_grid_matches_scalar("robinson", seed = 6613L)
})

test_that("cate_dnn grid diagnostics use fixed top-5 shares on every pseudo-outcome path", {
  tun <- .layer6_cate_tuning(default_s = 5L, s_grid = c(2L, 5L, 5L, 7L))
  seeds <- c(direct = 6619L, aipw = 6620L, robinson = 6621L)
  for (kind in names(seeds)) {
    case <- .layer6_cate_case(kind = kind, n = 40L, m = 3L, seed = seeds[[kind]])
    grid <- .layer6_do_cate(
      case$grid_args,
      method = "exact", k = 2L, tuning = tun, tuning_mode = "grid",
      se = FALSE,
      return = character(0)
    )

    .layer6_expect_fixed_top5_diagnostics(grid)
  }
})

test_that("cate_dnn grid mode keeps nuisance extras anchored to the default candidate", {
  tun <- .layer6_cate_tuning()
  for (kind in c("aipw", "robinson")) {
    case <- .layer6_cate_case(kind, seed = if (identical(kind, "aipw")) 6617L else 6618L)
    scalar <- .layer6_do_cate(
      case$scalar_args,
      method = "exact", k = 1L, tuning = tun,
      se = FALSE,
      return = c("psi", "nuisance", "tuning")
    )
    grid <- .layer6_do_cate(
      case$grid_args,
      method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
      se = FALSE,
      return = c("psi", "nuisance", "tuning")
    )

    expect_equal(grid$psi, scalar$psi, tolerance = 0)
    expect_equal(grid$nuisance, scalar$nuisance, tolerance = 0)
    expect_equal(grid$meta$tuning, scalar$meta$tuning)
  }
})

test_that("cate_dnn grid mode rejects explicit s and can use the paper default grid", {
  case <- .layer6_cate_case("direct", n = 40L, m = 2L, seed = 6614L)

  expect_error(
    .layer6_do_cate(
      case$grid_args,
      method = "exact", k = 1L, s = 4L,
      tuning_mode = "grid", se = FALSE
    ),
    "tuning_mode.*grid.*s"
  )

  grid <- .layer6_do_cate(
    case$grid_args,
    method = "exact", k = 1L,
    tuning_mode = "grid", se = FALSE,
    return = character(0)
  )

  expect_true(is.list(grid$tuning_grid))
  expect_equal(grid$tuning_grid$candidates$lambda, c(0.75, 1, 1.25, 1.5))
  expect_equal(grid$tuning_grid$meta$name, "paper_undersmooth")
  expect_equal(grid$meta$tuning$source, "paper_default")
})

test_that("cate_dnn grid mode constructs nuisance pseudo-outcomes once per fit", {
  aipw_case <- .layer6_cate_case("aipw", n = 40L, m = 2L, seed = 6615L)
  robinson_case <- .layer6_cate_case("robinson", n = 40L, m = 2L, seed = 6616L)
  tun <- .layer6_cate_tuning()

  aipw_calls <- .layer6_count_traced_calls("cate_rho_aipw", {
    grid <- .layer6_do_cate(
      aipw_case$grid_args,
      method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
      se = FALSE,
      return = "psi"
    )
    expect_true(is.list(grid$tuning_grid))
    expect_true(is.numeric(grid$psi))
  })
  expect_equal(aipw_calls, 1L)

  robinson_calls <- .layer6_count_traced_calls("cate_rho_rlearner", {
    grid <- .layer6_do_cate(
      robinson_case$grid_args,
      method = "exact", k = 1L, tuning = tun, tuning_mode = "grid",
      se = FALSE,
      return = "psi"
    )
    expect_true(is.list(grid$tuning_grid))
    expect_true(is.numeric(grid$psi))
  })
  expect_equal(robinson_calls, 1L)
})

test_that("later layers add tuning_mode to all main DNN wrappers", {
  expect_true("tuning_mode" %in% names(formals(np_dnn)))
  expect_true("tuning_mode" %in% names(formals(cate_dnn)))
  expect_true("tuning_mode" %in% names(formals(plqr_dnn)))
})
