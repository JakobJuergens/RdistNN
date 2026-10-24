.grid_plot_tuning <- function(default_s = 5L, s_grid = c(1L, 3L, 5L, 7L)) {
  new_tuning(
    name = "layer10_grid",
    eval = function(...) list(s = default_s, meta = list(lambda = 1)),
    eval_grid = function(...) {
      data.frame(
        lambda = c(0.4, 0.7, 1, 1.4),
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
}

.grid_plot_np_data <- function(n = 46L, m = 5L, seed = 1010L) {
  .test_rng_scope({
    X <- matrix(rnorm(n * 2L), n, 2L,
                dimnames = list(NULL, c("X1", "X2")))
    X0 <- cbind(seq(-1, 1, length.out = m), 0)
    y <- X[, 1L] - 0.25 * X[, 2L] + rnorm(n, sd = 0.2)
    list(X = X, X0 = X0, y = y)
  }, seed = seed)
}

.grid_plot_np_grid <- function(return = c("W", "ranks", "dx"),
                             se = TRUE, seed = 1010L) {
  d <- .grid_plot_np_data(seed = seed)
  .with_test_warning_families(np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 2L, tuning = .grid_plot_tuning(),
    tuning_mode = "grid", se = se, jackknife = "fast",
    return = return
  ), "infeasible_candidates")
}

.grid_plot_cate_grid <- function(seed = 1011L) {
  d <- .grid_plot_np_data(seed = seed)
  psi <- d$y - mean(d$y)
  .with_test_warning_families(cate_dnn(
    X = d$X, X0 = d$X0, pseudo_outcome = psi,
    method = "exact", k = 2L, tuning = .grid_plot_tuning(),
    tuning_mode = "grid", se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "dx", "psi")
  ), "infeasible_candidates")
}

.grid_plot_plqr_grid <- function(seed = 1012L) {
  demo <- plqr_demo_data(n = 50L, p = 2L, p_d = 2L, tau = 0.5, seed = seed)
  df <- demo$data
  .with_test_warning_families(plqr_dnn(
    X = as.matrix(df[paste0("X", 1:2)]),
    D = as.matrix(df[paste0("D", 1:2)]),
    g = as.numeric(df$g_oracle),
    mD = as.matrix(df[paste0("mD", 1:2)]),
    Y = as.numeric(df$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = cbind(seq(-0.5, 0.5, length.out = 4L), 0),
    method = "exact", k = 2L, tuning = .grid_plot_tuning(),
    tuning_mode = "grid", se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "dx")
  ), "infeasible_candidates")
}

test_that("grid fixtures filter intended warning families and preserve unrelated warnings", {
  for (case in list(
    list(fixture = .grid_plot_np_grid, estimator = "np_dnn"),
    list(fixture = .grid_plot_cate_grid, estimator = "cate_dnn"),
    list(fixture = .grid_plot_plqr_grid, estimator = "plqr_dnn")
  )) {
    real_estimator <- get(case$estimator, mode = "function")
    sentinel <- simpleWarning("unrelated fixture warning")
    fixture <- case$fixture
    environment(fixture) <- list2env(stats::setNames(list(function(...) {
      warning(sentinel)
      real_estimator(...)
    }), case$estimator), parent = environment(fixture))
    seen <- list()
    withCallingHandlers(fixture(), warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    expect_equal(length(seen), 1L, info = case$estimator)
    if (length(seen)) expect_identical(seen[[1L]], sentinel, info = case$estimator)
  }
})

test_that("diagnostics exposes candidate-level grid diagnostics on request", {
  skip_if_not_installed("ggplot2")
  fit <- .grid_plot_np_grid()

  d_default <- diagnostics(fit)
  expect_false("candidate_id" %in% names(d_default))
  expect_equal(nrow(d_default), nrow(fit$x0))

  d_all <- diagnostics(fit, candidates = "all")
  expect_s3_class(d_all, "rdistnn_diagnostics")
  expect_true(all(c(
    "candidate_id", "candidate_label", "is_default", "status", "reason",
    "n_eff", "max_weight", "top5_weight_share"
  ) %in% names(d_all)))
  expect_equal(d_all$candidate_id, fit$tuning_grid$candidates$candidate_id)
  expect_equal(d_all$status, fit$tuning_grid$candidates$status)
  expect_true(any(d_all$status != "ok"))
  expect_true(all(is.na(d_all$n_eff[d_all$status != "ok"])))

  ok_idx <- which(d_all$status == "ok")[[1L]]
  cand <- fit$tuning_grid$candidates[
    fit$tuning_grid$candidates$candidate_id == d_all$candidate_id[[ok_idx]],
    ,
    drop = FALSE
  ]
  tab <- rank_weight_dnn_exact(k = fit$meta$k, s = cand$s[[1L]])$table(fit$meta$n)
  denom <- sum(abs(tab))
  expected_top5 <- sum(abs(tab[seq_len(min(5L, length(tab)))])) / denom
  estimator_k_share <- sum(abs(tab[seq_len(min(fit$meta$k, length(tab)))])) / denom
  expect_equal(d_all$top5_weight_share[[ok_idx]], expected_top5, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(expected_top5, estimator_k_share)))

  p <- ggplot2::autoplot(d_all)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$candidate_label),
               d_all$candidate_label[d_all$status == "ok"])
  expect_match(p$labels$caption, "Dropped 1 unavailable candidate")
})

test_that("candidate estimate plots use tuning_grid rows, labels, and compact drops", {
  skip_if_not_installed("ggplot2")
  fit <- .grid_plot_np_grid()

  p <- plot_estimates(fit, candidates = "all")
  expect_s3_class(p, "ggplot")
  expect_true("candidate_id" %in% names(p$data))
  expect_identical(p$labels$colour, "candidate")
  expect_identical(p$labels$fill, "candidate")
  expect_match(p$labels$caption, "Dropped 1 unavailable candidate")
  expect_true("GeomLine" %in% .plot_layer_geoms(p))
  expect_equal(
    levels(p$data$candidate_label),
    c("lambda = 0.7 (s = 3)", "lambda = 1 (s = 5)", "lambda = 1.4 (s = 7)")
  )
  expect_identical(tail(unique(p$data$candidate_id), 1L), "lambda_1")

  scalar <- np_dnn(
    X = .grid_plot_np_data()$X,
    X0 = .grid_plot_np_data()$X0,
    Y = .grid_plot_np_data()$y,
    method = "exact", k = 2L, s = 5L, se = FALSE
  )
  expect_error(
    plot_estimates(scalar, candidates = "all"),
    "tuning_mode = \"grid\""
  )
  expect_error(
    plot_estimates(fit, candidates = "not_a_candidate"),
    "Available candidate ids"
  )
})

test_that("candidate weight plots compare rank profiles without stored W", {
  skip_if_not_installed("ggplot2")
  fit <- .grid_plot_np_grid(return = "ranks", se = FALSE)

  p <- plot_weights(fit, type = "profile", candidates = "all")
  expect_s3_class(p, "ggplot")
  expect_true("candidate_id" %in% names(p$data))
  expect_identical(p$labels$colour, "candidate")
  expect_identical(p$labels$linetype, "candidate")
  expect_match(p$labels$caption, "Dropped 1 unavailable candidate")

  expect_error(
    plot_weights(fit, type = "kernel", candidates = "all"),
    "candidate comparisons"
  )
})

test_that("candidate weight plots materialize weights only for displayed targets", {
  skip_if_not_installed("ggplot2")
  untied <- .grid_plot_np_grid(return = "ranks", se = FALSE)
  X <- matrix(c(-3, -2, -1, 1, 2, 3, 4, 5))
  X0 <- matrix(c(0, .2, 1))
  tuning <- new_tuning("display_groups", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 5L)))
  tied <- .with_distance_tie_warning(np_dnn(X = X, X0 = X0, Y = seq_len(8L),
    tuning = tuning, tuning_mode = "grid", se = FALSE, return = "ranks"))
  real_average <- RdistNN:::.dnn_average_distance_weights
  weighted_columns <- integer()
  local_mocked_bindings(.dnn_average_distance_weights = function(W, ...) {
    weighted_columns <<- c(weighted_columns, ncol(W))
    real_average(W, ...)
  }, .dnn_compute_group_facts = function(...) stop("Prepared group facts must be reused."),
    .package = "RdistNN")

  single <- plot_weights(untied, candidates = "all")
  expect_true(all(weighted_columns == 1L))
  expect_length(weighted_columns, sum(untied$tuning_grid$candidates$status == "ok"))
  expect_equal(nrow(single$data), untied$meta$n * sum(untied$tuning_grid$candidates$status == "ok"))
  weighted_columns <- integer()
  selected <- .with_distance_tie_warning(plot_weights(tied, candidates = "all", targets = c(3L, 1L)))
  expect_true(all(weighted_columns == 2L))
  expect_length(weighted_columns, 3L)
  expect_setequal(selected$data$target_id, c(3L, 1L))
  for (s in unique(selected$data$s)) {
    expected <- .block_oracle_weights(X, X0, .block_exact_table(8L, 1L, s))
    for (j in c(3L, 1L)) {
      rows <- selected$data[selected$data$s == s & selected$data$target_id == j, ]
      rows <- rows[order(rows$rank), ]
      ord <- order(abs(X[, 1L] - X0[j, 1L]), seq_len(nrow(X)))
      expect_equal(rows$weight, expected[ord, j], tolerance = 1e-14)
    }
  }
})

test_that("plot_tuning uses stored grids or internal engines without public refits", {
  skip_if_not_installed("ggplot2")
  d <- .grid_plot_np_data(seed = 1020L)
  grid <- .grid_plot_np_grid(seed = 1020L)
  scalar <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 2L, tuning = .grid_plot_tuning(),
    se = FALSE, return = character(0)
  )

  grid_calls <- .count_package_calls("np_dnn", {
    p_grid <- plot_tuning(grid, se = FALSE)
  })
  expect_identical(grid_calls, 0L)
  expect_s3_class(p_grid, "ggplot")
  expect_true(all(c("candidate_id", "candidate_label", "selected") %in% names(p_grid$data)))
  expect_match(p_grid$labels$caption, "Dropped 1 unavailable candidate")

  scalar_calls <- .count_package_calls("np_dnn", {
    expect_warning(p_scalar <- plot_tuning(scalar, se = FALSE),
      class = "rdistnn_infeasible_candidates_warning")
  })
  expect_identical(scalar_calls, 0L)
  expect_equal(p_scalar$data$estimate, p_grid$data$estimate)

  cate_calls <- .count_package_calls("cate_dnn", {
    p_cate <- plot_tuning(.grid_plot_cate_grid(), se = FALSE)
  })
  expect_identical(cate_calls, 0L)
  expect_true("GeomPointrange" %in% .plot_layer_geoms(p_cate))

  p_plqr <- plot_tuning(.grid_plot_plqr_grid(), se = FALSE)
  expect_s3_class(p_plqr, "ggplot")
  expect_identical(class(p_plqr$facet)[1L], "FacetWrap")
})

test_that("scalar PLQR tuning plots apply one requested uncertainty policy", {
  .with_distance_tie_warning({
  skip_if_not_installed("ggplot2")
  tuning <- new_tuning("scalar_plqr_plot_policy",
    eval = function(...) list(s = 2L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 4L)))
  args <- list(X = matrix(1:5), D = matrix(c(1, 2, 3, 2, 1)),
    mD = matrix(0, 5L, 1L), Y = c(-1, 1, -2, 2, 0), g = rep(0, 5L),
    theta_eval = .2, x0 = 2, tuning = tuning)
  for (fitted_se in c(FALSE, TRUE)) {
    fit <- do.call(plqr_dnn, c(args, list(se = fitted_se)))
    with_se <- plot_tuning(fit, se = TRUE)$data
    without_se <- plot_tuning(fit, se = FALSE)$data
    expect_true(all(is.finite(with_se$se)))
    expect_true(all(is.na(without_se$se)))
    expect_equal(with_se$estimate, without_se$estimate, tolerance = 1e-12)
    # A supplied fit-level SE is not a cache for a different plot policy.
    for (i in seq_len(nrow(with_se))) {
      scalar_args <- args
      scalar_args$tuning <- NULL
      scalar_args$s <- with_se$s[[i]]
      scalar <- do.call(plqr_dnn, scalar_args)
      expect_equal(with_se$se[[i]], unname(scalar$moment_se[1L, 1L]),
        tolerance = 1e-12)
    }

    stored <- do.call(plqr_dnn, c(args, list(se = fitted_se, tuning_mode = "grid")))
    stored_true <- plot_tuning(stored, se = TRUE)$data
    stored_false <- plot_tuning(stored, se = FALSE)$data
    expect_identical(stored_true$se, stored_false$se)
    expect_identical(is.finite(stored_true$se), rep(fitted_se, nrow(stored_true)))
  }
  })
})

test_that("scalar tuning feasibility follows the requested grid uncertainty", {
  skip_if_not_installed("ggplot2")
  tuning <- new_tuning("plot_order_boundary", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 8L)))
  for (constructor in list(function(...) np_dnn(Y = (1:8)^2, ...),
      function(...) cate_dnn(pseudo_outcome = (1:8)^2, ...))) {
    fit <- constructor(X = matrix(1:8), x0 = .25, tuning = tuning, se = TRUE)
    point <- plot_tuning(fit, se = FALSE)$data
    expect_true(8L %in% point$s)
    expect_true(all(is.na(point$se)))
    expect_identical(fit$meta$jackknife, "fast")
    expect_true(all(is.finite(fit$se)))
  }
})

test_that("grid dashboards include a tuning panel while preserving default panels", {
  skip_if_not_installed("ggplot2")
  fit <- .grid_plot_np_grid()

  testthat::local_mocked_bindings(
    .plot_pkg_installed = function(pkg) pkg != "patchwork",
    .package = "RdistNN"
  )
  suppressMessages(panels <- plot_dashboard(fit, level = 0.8))
  expect_true("tuning" %in% names(panels))
  expect_identical(attr(panels$estimates$data, "level"), 0.8)
  expect_false("candidate_id" %in% names(panels$estimates$data))
})
