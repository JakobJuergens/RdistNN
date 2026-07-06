# Helper-level tests: public wrappers do not call the grid engine in layer 4.
skip_if_not(
  all(vapply(
    c(
      "new_tuning",
      "rank_weight_dnn_exact",
      "rank_weight_eval_batch",
      "as_rank_weight"
    ),
    exists, logical(1),
    mode = "function"
  )),
  "DNN grid-engine dependencies not found; skipping file"
)

.layer4_ranks <- function() {
  R <- matrix(
    as.integer(c(
      1L, 2L, 3L, 4L, 5L,
      5L, 4L, 3L, 2L, 1L,
      2L, 3L, 4L, 5L, 1L
    )),
    nrow = 5L,
    ncol = 3L
  )
  R
}

.layer4_candidates <- function(
  s_grid = c(1L, 2L, 2L, 4L),
  kappa_grid = c(0.5, 1, 1.25, 2),
  default_s = 2L,
  default_kappa = 1,
  n = 5L,
  k = 2L,
  se = TRUE,
  jackknife = "fast"
) {
  tun <- new_tuning(
    name = "layer4_probe",
    eval = function(...) list(s = default_s, meta = list(kappa = default_kappa)),
    eval_grid = function(...) {
      data.frame(
        kappa = kappa_grid,
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
  RdistNN:::.dnn_resolve_s(
    s = NULL,
    tuning = tun,
    n = n,
    p = 2L,
    k = k,
    se = se,
    jackknife = jackknife
  )$tuning$grid$candidates
}

.layer4_estimator <- function(signal) {
  force(signal)
  function(W, candidate_row) {
    data.frame(
      target_id = seq_len(ncol(W)),
      component = NA_character_,
      estimate = as.numeric(crossprod(signal, W)),
      stringsAsFactors = FALSE
    )
  }
}

.layer4_se <- function(ranks) {
  force(ranks)
  function(candidate_row, tab_jk) {
    data.frame(
      target_id = seq_len(ncol(ranks)),
      component = NA_character_,
      se = rep(sum(tab_jk), ncol(ranks)),
      stringsAsFactors = FALSE
    )
  }
}

.layer4_counted_dnn_spec <- function(calls, fail_s = integer(0L)) {
  force(calls)
  force(fail_s)
  function(candidate_row) {
    s <- as.integer(candidate_row$s[[1L]])
    calls$keys <- c(calls$keys, as.character(candidate_row$computed_key[[1L]]))
    if (s %in% fail_s) {
      stop(sprintf("boom s=%d", s), call. = FALSE)
    }
    rank_weight_dnn_exact(k = 1L, s = s)
  }
}

.layer4_fixed_table_spec <- function(tab) {
  force(tab)
  spec <- list(
    name = "fixed_table",
    table = function(n) {
      if (n == length(tab)) {
        return(tab)
      }
      rep(1 / n, n)
    },
    .__kind = "rank_weight"
  )
  class(spec) <- c("rank_weight", class(spec))
  spec
}

test_that("grid engine computes duplicate feasible keys once and reports all rows", {
  ranks <- .layer4_ranks()
  signal <- c(2, -1, 3, 0, 4)
  candidates <- .layer4_candidates()
  calls <- new.env(parent = emptyenv())
  calls$keys <- character()

  out <- RdistNN:::.dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = .layer4_counted_dnn_spec(calls),
    estimate_fn = .layer4_estimator(signal),
    se_fn = .layer4_se(ranks),
    ctx = list(k = 2L)
  )

  expect_named(out, c("candidates", "estimates", "diagnostics", "meta"))
  expect_equal(calls$keys, c("s_2", "s_4"))
  expect_equal(out$candidates$candidate_id, candidates$candidate_id)
  expect_equal(out$candidates$status, c("infeasible", "ok", "ok", "ok"))
  expect_equal(nrow(out$estimates), nrow(candidates) * ncol(ranks))
  expect_false("W" %in% names(out))

  dup <- out$estimates[out$estimates$computed_key == "s_2", ]
  expect_equal(unique(dup$status), "ok")
  expect_equal(
    dup$estimate[dup$candidate_id == "kappa_1"],
    dup$estimate[dup$candidate_id == "kappa_1p25"],
    tolerance = 1e-12
  )
  expect_equal(out$meta$n_feasible, 3L)
  expect_equal(out$meta$n_skipped, 1L)
  expect_equal(out$meta$n_unique, 2L)
})

test_that("grid engine emits NA rows for skipped infeasible candidates", {
  ranks <- .layer4_ranks()
  candidates <- .layer4_candidates()
  out <- RdistNN:::.dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = .layer4_counted_dnn_spec(new.env(parent = emptyenv())),
    estimate_fn = .layer4_estimator(c(1, 2, 3, 4, 5)),
    se_fn = .layer4_se(ranks),
    ctx = list(k = 2L)
  )

  skipped_id <- candidates$candidate_id[!candidates$feasible]
  skipped_est <- out$estimates[out$estimates$candidate_id == skipped_id, ]
  skipped_diag <- out$diagnostics[out$diagnostics$candidate_id == skipped_id, ]

  expect_equal(skipped_est$status, rep("infeasible", ncol(ranks)))
  expect_true(all(is.na(skipped_est$estimate)))
  expect_true(all(is.na(skipped_est$se)))
  expect_equal(skipped_diag$status, "infeasible")
  expect_true(is.na(skipped_diag$ess))
  expect_equal(skipped_diag$top_n, 5L)
  expect_true(is.na(skipped_diag$top5_weight_share))
  expect_true(is.na(skipped_diag$top_share))
})

test_that("grid engine records failed non-default candidate rows", {
  ranks <- .layer4_ranks()
  candidates <- .layer4_candidates()
  calls <- new.env(parent = emptyenv())
  calls$keys <- character()

  out <- RdistNN:::.dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = .layer4_counted_dnn_spec(calls, fail_s = 4L),
    estimate_fn = .layer4_estimator(c(1, 2, 3, 4, 5)),
    se_fn = .layer4_se(ranks),
    ctx = list(k = 2L)
  )

  failed_id <- candidates$candidate_id[candidates$s == 4L]
  failed_est <- out$estimates[out$estimates$candidate_id == failed_id, ]
  expect_equal(out$candidates$status[out$candidates$candidate_id == failed_id], "failed")
  expect_match(out$candidates$reason[out$candidates$candidate_id == failed_id], "boom s=4")
  expect_equal(failed_est$status, rep("failed", ncol(ranks)))
  expect_true(all(is.na(failed_est$estimate)))
  expect_true(all(is.na(failed_est$se)))
  expect_equal(out$meta$n_failed, 1L)
})

test_that("grid engine rethrows failed default candidates", {
  ranks <- .layer4_ranks()
  candidates <- .layer4_candidates()
  expect_error(
    RdistNN:::.dnn_grid_eval(
      candidates = candidates,
      ranks = ranks,
      weight_spec_fn = .layer4_counted_dnn_spec(new.env(parent = emptyenv()), fail_s = 2L),
      estimate_fn = .layer4_estimator(c(1, 2, 3, 4, 5)),
      se_fn = .layer4_se(ranks),
      ctx = list(k = 2L)
    ),
    "boom s=2"
  )
})

test_that("grid engine computes table diagnostics from the rank-weight table", {
  ranks <- .layer4_ranks()
  candidates <- .layer4_candidates(
    s_grid = 3L,
    kappa_grid = 1,
    default_s = 3L,
    n = 5L,
    k = 2L,
    se = FALSE,
    jackknife = "none"
  )
  tab <- c(0.5, 0.3, 0.2, 0, 0)

  out <- RdistNN:::.dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = function(candidate_row) .layer4_fixed_table_spec(tab),
    estimate_fn = .layer4_estimator(c(1, 2, 3, 4, 5)),
    se_fn = NULL,
    ctx = list(k = 2L)
  )

  expect_equal(out$diagnostics$ess, 1 / sum(tab^2), tolerance = 1e-12)
  expect_equal(out$diagnostics$max_weight, max(abs(tab)), tolerance = 1e-12)
  expect_equal(out$diagnostics$top_n, 5L)
  expect_equal(out$diagnostics$top5_weight_share, 1, tolerance = 1e-12)
  expect_equal(out$diagnostics$top_share, out$diagnostics$top5_weight_share, tolerance = 1e-12)
  expect_false(out$diagnostics$degenerate)
})

test_that("grid engine diagnostic top share is independent of estimator k", {
  ranks <- matrix(as.integer(seq_len(6L)), ncol = 1L)
  candidates <- .layer4_candidates(
    s_grid = 3L,
    kappa_grid = 1,
    default_s = 3L,
    n = 6L,
    k = 1L,
    se = FALSE,
    jackknife = "none"
  )
  tab <- c(0.4, 0.2, 0.1, 0.1, 0.1, 0.1)

  out <- RdistNN:::.dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = function(candidate_row) .layer4_fixed_table_spec(tab),
    estimate_fn = .layer4_estimator(seq_len(6L)),
    se_fn = NULL,
    ctx = list(k = 1L)
  )

  expect_equal(out$diagnostics$top_n, 5L)
  expect_equal(out$diagnostics$top5_weight_share, 0.9, tolerance = 1e-12)
  expect_equal(out$diagnostics$top_share, out$diagnostics$top5_weight_share, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(out$diagnostics$top5_weight_share, tab[[1L]])))
})

test_that("grid engine falls back to rank_weight_eval_batch when no table is available", {
  ranks <- .layer4_ranks()
  candidates <- .layer4_candidates(
    s_grid = 2L,
    kappa_grid = 1,
    default_s = 2L,
    n = 5L,
    k = 1L,
    se = FALSE,
    jackknife = "none"
  )
  signal <- c(1, -1, 2, -2, 3)
  spec <- as_rank_weight(function(r) {
    w <- 6 - as.numeric(r)
    w / sum(w)
  }, name = "layer4_no_table")
  expected_W <- rank_weight_eval_batch(spec, ranks)

  out <- RdistNN:::.dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = function(candidate_row) spec,
    estimate_fn = .layer4_estimator(signal),
    se_fn = NULL,
    ctx = list(k = 2L)
  )

  expect_equal(out$estimates$estimate, as.numeric(crossprod(signal, expected_W)), tolerance = 1e-12)
  expect_true(all(is.na(out$diagnostics$ess)))
  expect_true(all(is.na(out$diagnostics$top5_weight_share)))
  expect_true(all(is.na(out$diagnostics$top_share)))
})

test_that("later layers wire plqr wrappers while plotting stays outside the grid engine", {
  expect_true("tuning_mode" %in% names(formals(plqr_dnn)))

  source_dir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
  skip_if_not(dir.exists(source_dir), "source R directory not available in installed-package checks")
  source_refs <- unlist(lapply(
    list.files(source_dir, pattern = "[.]R$", full.names = TRUE),
    function(path) {
      txt <- readLines(path, warn = FALSE)
      paste(path, grep("[.]dnn_grid_eval", txt, value = TRUE), sep = ":")
    }
  ), use.names = FALSE)
  source_refs <- source_refs[grepl("[.]dnn_grid_eval", source_refs)]
  expect_true(any(grepl("user-dnn-core[.]R", source_refs)))
  expect_true(any(grepl("user-dnn-plqr[.]R", source_refs)))
  expect_false(any(grepl("plot-|sim-", source_refs)))
})
