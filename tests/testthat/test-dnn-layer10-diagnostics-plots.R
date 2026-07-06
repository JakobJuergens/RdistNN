skip_if_not(
  all(vapply(
    c("np_dnn", "cate_dnn", "plqr_dnn", "new_tuning", "diagnostics",
      "plot_estimates", "plot_weights", "plot_tuning", "plot_dashboard"),
    exists, logical(1),
    mode = "function"
  )),
  "Layer 10 diagnostics/plot dependencies not found; skipping file"
)

.layer10_geoms <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
}

.layer10_count_traced_calls <- function(fun, expr) {
  trace_env <- new.env(parent = emptyenv())
  trace_env$calls <- 0L
  assign(".layer10_trace_env", trace_env, envir = .GlobalEnv)
  on.exit(rm(".layer10_trace_env", envir = .GlobalEnv), add = TRUE)
  suppressMessages(invisible(trace(
    fun,
    tracer = quote({
      .GlobalEnv$.layer10_trace_env$calls <- .GlobalEnv$.layer10_trace_env$calls + 1L
    }),
    where = asNamespace("RdistNN"),
    print = FALSE
  )))
  on.exit(suppressMessages(untrace(fun, where = asNamespace("RdistNN"))), add = TRUE)
  force(expr)
  trace_env$calls
}

.layer10_tuning <- function(default_s = 5L, s_grid = c(1L, 3L, 5L, 7L)) {
  new_tuning(
    name = "layer10_grid",
    eval = function(...) list(s = default_s, meta = list(kappa = 1)),
    eval_grid = function(...) {
      data.frame(
        kappa = c(0.4, 0.7, 1, 1.4),
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
}

.layer10_np_data <- function(n = 46L, m = 5L, seed = 1010L) {
  set.seed(seed)
  X <- matrix(rnorm(n * 2L), n, 2L,
              dimnames = list(NULL, c("X1", "X2")))
  X0 <- cbind(seq(-1, 1, length.out = m), 0)
  y <- X[, 1L] - 0.25 * X[, 2L] + rnorm(n, sd = 0.2)
  list(X = X, X0 = X0, y = y)
}

.layer10_np_grid <- function(return = c("W", "ranks", "dx"),
                             se = TRUE, seed = 1010L) {
  d <- .layer10_np_data(seed = seed)
  suppressWarnings(np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 2L, tuning = .layer10_tuning(),
    tuning_mode = "grid", se = se, jackknife = "fast",
    return = return
  ))
}

.layer10_cate_grid <- function(seed = 1011L) {
  d <- .layer10_np_data(seed = seed)
  psi <- d$y - mean(d$y)
  suppressWarnings(cate_dnn(
    X = d$X, X0 = d$X0, pseudo_outcome = psi,
    method = "exact", k = 2L, tuning = .layer10_tuning(),
    tuning_mode = "grid", se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "dx", "psi")
  ))
}

.layer10_plqr_grid <- function(seed = 1012L) {
  demo <- plqr_demo_data(n = 50L, p = 2L, p_d = 2L, tau = 0.5, seed = seed)
  df <- demo$data
  suppressWarnings(plqr_dnn(
    X = as.matrix(df[paste0("X", 1:2)]),
    D = as.matrix(df[paste0("D", 1:2)]),
    q = as.numeric(df$q_oracle),
    mD = as.matrix(df[paste0("mD", 1:2)]),
    Y = as.numeric(df$Y),
    beta = demo$truth$beta_true, tau = 0.5,
    X0 = cbind(seq(-0.5, 0.5, length.out = 4L), 0),
    method = "exact", k = 2L, tuning = .layer10_tuning(),
    tuning_mode = "grid", se = TRUE, jackknife = "fast",
    return = c("W", "ranks", "dx")
  ))
}

test_that("diagnostics exposes candidate-level grid diagnostics on request", {
  skip_if_not_installed("ggplot2")
  fit <- .layer10_np_grid()

  d_default <- diagnostics(fit)
  expect_false("candidate_id" %in% names(d_default))
  expect_equal(nrow(d_default), nrow(fit$x0))

  d_all <- diagnostics(fit, candidates = "all")
  expect_s3_class(d_all, "dnn_diagnostics")
  expect_true(all(c(
    "candidate_id", "candidate_label", "is_default", "status", "reason",
    "n_eff", "max_weight", "top5_weight_share"
  ) %in% names(d_all)))
  expect_equal(d_all$candidate_id, fit$tuning_grid$candidates$candidate_id)
  expect_equal(d_all$status, fit$tuning_grid$candidates$status)
  expect_true(any(d_all$status != "ok"))
  expect_true(all(is.na(d_all$n_eff[d_all$status != "ok"])))

  p <- ggplot2::autoplot(d_all)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$candidate_label),
               d_all$candidate_label[d_all$status == "ok"])
  expect_match(p$labels$caption, "Dropped 1 unavailable candidate")
})

test_that("candidate estimate plots use tuning_grid rows, labels, and compact drops", {
  skip_if_not_installed("ggplot2")
  fit <- .layer10_np_grid()

  p <- plot_estimates(fit, candidates = "all")
  expect_s3_class(p, "ggplot")
  expect_true("candidate_id" %in% names(p$data))
  expect_identical(p$labels$colour, "candidate")
  expect_identical(p$labels$fill, "candidate")
  expect_match(p$labels$caption, "Dropped 1 unavailable candidate")
  expect_true("GeomLine" %in% .layer10_geoms(p))
  expect_equal(
    levels(p$data$candidate_label),
    c("kappa = 0.7 (s = 3)", "kappa = 1 (s = 5)", "kappa = 1.4 (s = 7)")
  )
  expect_identical(tail(unique(p$data$candidate_id), 1L), "kappa_1")

  scalar <- np_dnn(
    X = .layer10_np_data()$X,
    X0 = .layer10_np_data()$X0,
    Y = .layer10_np_data()$y,
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
  fit <- .layer10_np_grid(return = "ranks", se = FALSE)

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

test_that("plot_tuning uses stored grids or internal engines without public refits", {
  skip_if_not_installed("ggplot2")
  d <- .layer10_np_data(seed = 1020L)
  grid <- .layer10_np_grid(seed = 1020L)
  scalar <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "exact", k = 2L, tuning = .layer10_tuning(),
    se = FALSE, return = character(0)
  )

  grid_calls <- .layer10_count_traced_calls("np_dnn", {
    p_grid <- plot_tuning(grid, se = FALSE)
  })
  expect_identical(grid_calls, 0L)
  expect_s3_class(p_grid, "ggplot")
  expect_true(all(c("candidate_id", "candidate_label", "selected") %in% names(p_grid$data)))
  expect_match(p_grid$labels$caption, "Dropped 1 unavailable candidate")

  scalar_calls <- .layer10_count_traced_calls("np_dnn", {
    p_scalar <- plot_tuning(scalar, se = FALSE)
  })
  expect_identical(scalar_calls, 0L)
  expect_equal(p_scalar$data$estimate, p_grid$data$estimate)

  cate_calls <- .layer10_count_traced_calls("cate_dnn", {
    p_cate <- plot_tuning(.layer10_cate_grid(), se = FALSE)
  })
  expect_identical(cate_calls, 0L)
  expect_true("GeomPointrange" %in% .layer10_geoms(p_cate))

  p_plqr <- plot_tuning(.layer10_plqr_grid(), se = FALSE)
  expect_s3_class(p_plqr, "ggplot")
  expect_identical(class(p_plqr$facet)[1L], "FacetWrap")
})

test_that("grid dashboards include a tuning panel while preserving default panels", {
  skip_if_not_installed("ggplot2")
  fit <- .layer10_np_grid()

  testthat::local_mocked_bindings(
    .plot_pkg_installed = function(pkg) pkg != "patchwork",
    .package = "RdistNN"
  )
  suppressMessages(panels <- plot_dashboard(fit, level = 0.8))
  expect_true("tuning" %in% names(panels))
  expect_identical(attr(panels$estimates$data, "level"), 0.8)
  expect_false("candidate_id" %in% names(panels$estimates$data))
})
