.retained_report_fits <- function(se = TRUE, method = "exact") {
  X <- matrix(seq_len(8L), ncol = 1L)
  X0 <- matrix(c(0, 9, 0), ncol = 1L)
  y <- seq_len(8L)^2
  tuning <- new_tuning("retained_report_state", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L)))
  common <- list(X = X, X0 = X0, se = se, method = method, tuning = tuning,
    tuning_mode = "grid", return = character())
  constructors <- list(
    np = function(...) np_dnn(Y = y, ...),
    cate = function(...) cate_dnn(pseudo_outcome = y, ...),
    ratio = function(...) affine_ratio_dnn(numerator = y, denominator = rep(1, 8L), ...),
    clate = function(...) clate_dnn(Y = y, D = rep(c(0, 1), 4L), A = rep(c(0, 1), 4L),
      e = rep(.5, 8L), m_y0 = rep(0, 8L), m_y1 = rep(1, 8L),
      m_d0 = rep(0, 8L), m_d1 = rep(1, 8L), ...))
  fits <- lapply(constructors, function(fn) .with_localization_policy_warnings(do.call(fn, common)))
  fits$plqr <- plqr_dnn(X = X, X0 = X0, Y = rep(1, 8L),
    D = matrix(seq_len(8L)), mD = matrix(0, 8L, 1L), g = rep(0, 8L),
    theta_eval = 0, density = rep(1, 8L), s = 3L, se = se, method = method, return = character())
  fits$tdnn <- np_tdnn(X = X, X0 = X0, Y = y, s1 = 2L, s2 = 3L,
    se = se, method = method, return = character())
  fits
}

.retained_report_routes <- function(fit) {
  routes <- list(coef = coef, summary = summary, confint = confint, vcov = vcov,
    estimates = plot_data_estimates, diagnostics = diagnostics,
    prediction = function(x) predict(x, x0 = .25))
  if (!is.null(fit$tuning_grid)) routes$candidate <- function(x) candidate_fit(x, 1L)
  routes
}

.retained_control_fits <- function(method = "exact") {
  fits <- .retained_report_fits(se = FALSE, method = method)
  tuning <- new_tuning("retained_control_state", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L)))
  fits$plqr <- plqr_dnn(X = matrix(seq_len(8L)), X0 = fits$np$x0,
    Y = rep(1, 8L), D = matrix(seq_len(8L)), mD = matrix(0, 8L),
    g = rep(0, 8L), theta_eval = 0, se = FALSE, method = method,
    tuning = tuning, tuning_mode = "grid", return = character())
  model <- new_gmm_model("retained_control_mean",
    setup = function(sample, n, ...) list(y = sample$y),
    moment_obs = function(sample, theta, ctx) ctx$y - theta[1L],
    start = function(sample, ctx, w, x0) sum(w * ctx$y), coef_names = "mean")
  fits$gmm <- gmm_dnn(X = matrix(seq_len(8L)), X0 = fits$np$x0,
    model_data = data.frame(y = seq_len(8L)^2), model = model,
    s = 3L, se = FALSE, method = method, return = character())
  fits
}

.retained_control_routes <- function(fit) {
  routes <- list(coef = coef)
  if (!inherits(fit, "rdistnn_gmm_result")) {
    routes$prediction <- function(x) predict(x, X0 = x$x0)
  }
  if (!is.null(fit$tuning_grid)) routes$candidate <- function(x) candidate_fit(x, 1L)
  routes
}

.retained_scalar_tuning_fits <- function(tuning = NULL, s = NULL) {
  X <- matrix(seq_len(12L))
  y <- seq_len(12L)^2
  common <- list(X = X, x0 = .25, tuning = tuning, s = s,
    se = FALSE, return = c("W", "tuning"))
  constructors <- list(
    np = function(...) np_dnn(Y = y, ...),
    cate = function(...) cate_dnn(pseudo_outcome = y, ...),
    ratio = function(...) affine_ratio_dnn(numerator = y, denominator = rep(1, 12L), ...))
  lapply(constructors, function(fn) do.call(fn, common))
}

test_that("scalar fits validate present current tuning metadata before reuse", {
  tuning <- new_tuning("retained_scalar_tuning",
    eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 4L)))
  fits <- c(.retained_scalar_tuning_fits(tuning), .retained_scalar_tuning_fits())
  changes <- list(
    function(x) { x$meta$tuning <- list(source = "old"); x },
    function(x) { x$meta$tuning$default$s <- 99L; x },
    function(x) { x$meta$tuning$default$meta <- NULL; x },
    function(x) { x$meta$tuning$grid <- NULL; x },
    function(x) { x$meta["tuning"] <- list(NULL); x })
  consumers <- list(coef, function(x) predict(x, x0 = .25), diagnostics)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    consumers <- c(consumers, list(function(x) plot_tuning(x, se = FALSE)))
  }
  for (fit in fits) for (change in changes) for (consume in consumers) {
    err <- tryCatch(consume(change(fit)), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    if (inherits(err, "condition")) {
      expect_match(conditionMessage(err), "tuning", fixed = TRUE)
      expect_match(conditionMessage(err), "Refit", fixed = TRUE)
    }
  }
})

test_that("current scalar tuning and manual optional records remain reusable", {
  tuning <- new_tuning("retained_scalar_passthrough",
    eval = function(...) list(s = 3L,
      meta = list(name = matrix("descriptive", 1L, 1L), n = 99L, kx = 7L)),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 4L)))
  fits <- c(.retained_scalar_tuning_fits(tuning), .retained_scalar_tuning_fits(),
    .retained_scalar_tuning_fits(s = 3L))
  for (fit in fits) {
    saved <- unserialize(serialize(fit, NULL))
    expect_equal(coef(saved), coef(fit))
    expect_equal(as.numeric(predict(saved, x0 = .25)), as.numeric(coef(fit)))
    expect_s3_class(diagnostics(saved), "rdistnn_diagnostics")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      expect_s3_class(plot_tuning(saved, se = FALSE), "ggplot")
    }
  }
  manual <- .retained_scalar_tuning_fits(s = 3L)$np
  expect_false("tuning" %in% names(manual$meta))
  expect_identical(diagnostics(manual)$tuning_source, "manual")
  manual$tuning$default$s <- 4L
  err <- tryCatch(diagnostics(manual), error = identity)
  expect_s3_class(err, "rdistnn_state_error")
  if (inherits(err, "condition")) {
    expect_match(conditionMessage(err), "object$tuning$default$s", fixed = TRUE)
  }
})

test_that("fresh structural distance specifications retain current geometry identity", {
  distance <- list(name = "structural_distance",
    eval = function(X, x0, ctx) abs(X[, 1L] - x0))
  X <- matrix(seq_len(8L))
  y <- seq_len(8L)^2
  common <- list(X = X, x0 = .25, distance = distance, se = FALSE)
  constructors <- list(
    function(...) np_dnn(Y = y, s = 3L, ...),
    function(...) cate_dnn(pseudo_outcome = y, s = 3L, ...),
    function(...) affine_ratio_dnn(numerator = y, denominator = rep(1, 8L), s = 3L, ...),
    function(...) np_tdnn(Y = y, s1 = 2L, s2 = 3L, ...))
  for (fn in constructors) {
    fit <- do.call(fn, common)
    for (field in c("distance", "inference_distance")) {
      expect_s3_class(fit$.model[[field]], "rdistnn_distance")
    }
    saved <- unserialize(serialize(fit, NULL))
    expect_equal(coef(saved), coef(fit))
    expect_equal(as.numeric(predict(saved, x0 = .25)), as.numeric(coef(fit)))
  }
  expect_null(attr(distance, "class"))
})

test_that("altered model orders cannot disagree with cached coefficients", {
  X <- matrix(seq_len(8L))
  y <- seq_len(8L)^2
  fits <- list(np_dnn(X = X, x0 = .25, Y = y, s = 3L, se = FALSE),
    cate_dnn(X = X, x0 = .25, pseudo_outcome = y, s = 3L, se = FALSE),
    affine_ratio_dnn(X = X, x0 = .25, numerator = y, denominator = rep(1, 8L), s = 3L, se = FALSE))
  # Independent k = 1, s = 3 exact subset probabilities at the ordered target.
  expected <- sum(c(21, 15, 10, 6, 3, 1, 0, 0) * y) / 56
  expect_identical(expected, 6.75)
  for (fit in fits) {
    expect_equal(as.numeric(coef(fit)), expected)
    expect_equal(as.numeric(predict(fit, x0 = .25)), expected)
    bad <- fit
    bad$.model$s <- 4L
    expect_error(coef(bad), class = "rdistnn_state_error")
    expect_error(predict(bad, x0 = .25), class = "rdistnn_state_error")
  }
})

test_that("current retained fields are exact names rather than partial aliases", {
  fits <- .retained_report_fits(se = FALSE)[c("np", "cate", "ratio")]
  for (fit in fits) {
    for (field in c("x_mat", "groups", "ranks", "dx", "tie_policy")) {
      bad <- fit
      names(bad$.model)[names(bad$.model) == field] <- paste0(field, "_old")
      for (route in list(stats::coef, function(x) predict(x, x0 = .25))) {
        err <- tryCatch(route(bad), error = identity)
        expect_s3_class(err, "rdistnn_state_error")
        if (inherits(err, "condition")) expect_match(conditionMessage(err), "Refit", fixed = TRUE)
      }
    }
    bad <- fit
    names(bad)[names(bad) == ".model"] <- ".model_old"
    expect_error(coef(bad), class = "rdistnn_state_error")
    expect_error(predict(bad, x0 = .25), class = "rdistnn_state_error")
  }
})

test_that("retained metric identities are current specifications", {
  for (fit in .retained_report_fits(se = FALSE)[c("np", "cate", "ratio")]) {
    for (old_class in list("dnn_distance", NULL)) {
      bad <- fit
      class(bad$.model$distance) <- old_class
      class(bad$.model$inference_distance) <- old_class
      expect_error(coef(bad), class = "rdistnn_state_error")
      expect_error(predict(bad, x0 = .25), class = "rdistnn_state_error")
    }
  }
})

test_that("duplicated fitted controls must agree across every current family", {
  fits <- .retained_control_fits()
  for (family in names(fits)) {
    fit <- fits[[family]]
    fields <- if (family == "tdnn") c("method", "s1", "s2", "dim") else c("method", "k", "s")
    if (family %in% c("ratio", "clate")) fields <- c(fields, "normalize")
    for (field in fields) for (location in c("meta", ".model")) {
      changed <- switch(field, method = "asymp", k = 2L, s = 4L,
        s1 = 1L, s2 = 4L, dim = 2L, normalize = FALSE)
      for (case in c("missing", "malformed", "mismatch")) {
        bad <- fit
        if (case == "missing") bad[[location]][field] <- NULL
        if (case == "malformed") bad[[location]][field] <- list(c("invalid", "control"))
        if (case == "mismatch") bad[[location]][field] <- list(changed)
        for (route in .retained_control_routes(bad)) {
          error <- tryCatch(route(bad), error = identity)
          expect_s3_class(error, "rdistnn_state_error")
          if (inherits(error, "rdistnn_state_error")) {
            expect_identical(error$code, "state_invalid")
            expect_match(conditionMessage(error), "[Rr]efit")
          }
        }
      }
    }
    for (value in list(NULL, "TRUE", c(TRUE, FALSE), NA)) {
      bad <- fit
      bad$.model["normalize"] <- list(value)
      expect_error(coef(bad), class = "rdistnn_state_error")
    }
  }
})

test_that("geometry labels preserve their distinct current meanings", {
  for (fit in .retained_control_fits()) {
    for (field in c("distance", "tie_rule")) {
      bad <- fit
      bad$meta[field] <- NULL
      expect_error(coef(bad), class = "rdistnn_state_error")
    }
    bad <- fit
    bad$meta$tie_rule <- "obsolete_policy"
    expect_error(coef(bad), class = "rdistnn_state_error")
    bad <- fit
    bad$.model$inference_distance$name <- "different_geometry"
    expect_error(coef(bad), class = "rdistnn_state_error")
  }
  alias <- "retained_control_geometry_alias"
  register_distance(alias, distance_l2())
  on.exit(rm(list = alias, envir = RdistNN:::.dist_registry), add = TRUE)
  fit <- np_dnn(X = matrix(seq_len(8L)), x0 = .25, Y = seq_len(8L)^2,
    distance = alias, s = 3L, se = FALSE)
  expect_identical(fit$meta$distance, alias)
  expect_identical(fit$.model$distance$name, "l2")
  register_distance(alias, distance_l1())
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(as.numeric(coef(restored)), 6.75)
  expect_equal(as.numeric(predict(restored, x0 = .25)), 6.75)
})

test_that("retained metric callback schemas fail at the shared state boundary", {
  fits <- .retained_report_fits()[c("np", "tdnn", "cate", "ratio", "clate")]
  for (fit in fits) for (field in c("distance", "inference_distance")) {
    for (change in c("missing", "invalid_eval", "invalid_setup", "prefix", "duplicate")) {
      bad <- fit
      spec <- bad$.model[[field]]
      if (change %in% c("missing", "prefix")) {
        evaluator <- spec$eval
        spec[c("eval", "eval_batch")] <- NULL
        if (change == "prefix") spec$eval_backup <- evaluator
      }
      if (change == "invalid_eval") spec$eval <- 1L
      if (change == "invalid_setup") spec$setup <- 1L
      if (change == "duplicate") {
        spec <- c(spec, list(eval = spec$eval))
        class(spec) <- class(bad$.model[[field]])
      }
      bad$.model[[field]] <- spec
      for (consume in list(coef, vcov, function(x) predict(x, x0 = .25))) {
        error <- tryCatch({ consume(bad); NULL }, error = identity)
        expect_s3_class(error, "rdistnn_state_error")
        if (inherits(error, "condition")) {
          expect_identical(error$code, "state_invalid")
          expect_match(conditionMessage(error), "[Rr]efit")
        }
      }
    }
  }
})

test_that("retained metric aliases cannot hide behind another valid evaluator", {
  calls <- 0L
  scalar <- function(X, x0, ctx) { calls <<- calls + 1L; abs(X[, 1L] - x0[1L]) }
  batch <- function(X, X0, ctx) { calls <<- calls + 1L; abs(outer(X[, 1L], X0[, 1L], "-")) }
  fit <- np_dnn(X = matrix(seq_len(8L)), x0 = .25, Y = seq_len(8L)^2,
    distance = new_distance("retained_aliases", eval = scalar, eval_batch = batch),
    s = 3L, jackknife = "loo")
  for (field in c("distance", "inference_distance")) {
    for (alias in c("eval_batch", "setup")) {
      bad <- fit
      spec <- bad$.model[[field]]
      if (alias == "setup" && !"setup" %in% names(spec)) spec["setup"] <- list(NULL)
      names(spec)[names(spec) == alias] <- paste0(alias, "_backup")
      alternate <- if (alias == "eval_batch") "eval" else "eval_batch"
      expect_true(is.function(spec[[alternate, exact = TRUE]]))
      expect_identical(class(spec), class(bad$.model[[field]]))
      bad$.model[[field]] <- spec
      calls <- 0L
      for (consume in list(coef, summary, vcov, function(x) predict(x, x0 = .25))) {
        error <- tryCatch({ consume(bad); NULL }, error = identity)
        expect_s3_class(error, "rdistnn_state_error")
        if (inherits(error, "condition")) {
          expect_identical(error$code, "state_invalid")
          expect_match(conditionMessage(error), "[Rr]efit")
        }
      }
      expect_identical(calls, 0L)
    }
  }
})

test_that("retained metric validation accepts current callbacks without evaluating them", {
  X <- matrix(seq_len(8L))
  X0 <- matrix(c(.25, 6.25))
  calls <- 0L
  scalar <- function(X, x0, ctx) { calls <<- calls + 1L; abs(X[, 1L] - x0[1L]) }
  batch <- function(X, X0, ctx) { calls <<- calls + 1L; abs(outer(X[, 1L], X0[, 1L], "-")) }
  for (spec in list(new_distance("retained_scalar", eval = scalar, setup = NULL),
      new_distance("retained_batch", eval_batch = batch, setup = NULL),
      list(name = "retained_scalar_structural", eval = scalar),
      list(name = "retained_batch_structural", eval_batch = batch))) {
    spec$custom_detail <- "retained extension metadata"
    spec$name_note <- "metadata beside exact name"
    if ("eval" %in% names(spec)) spec$eval_note <- "metadata beside exact eval"
    fit <- np_dnn(X = X, X0 = X0, Y = seq_len(8L)^2,
      distance = spec, ctx = NULL, s = 3L, jackknife = "loo")
    calls <- 0L
    expect_no_error(coef(fit))
    expect_no_error(summary(fit))
    expect_identical(calls, 0L)
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(coef(restored), coef(fit))
    expect_equal(vcov(restored), vcov(fit))
    expect_equal(predict(restored, X0 = X0, se.fit = TRUE),
      predict(fit, X0 = X0, se.fit = TRUE))
    expect_true("ctx" %in% names(fit$.model))
    expect_null(fit$.model$ctx)
  }
  fit <- np_dnn(X = X, X0 = X0, Y = seq_len(8L)^2, s = 3L, jackknife = "loo")
  sentinel <- simpleError("opaque retained metric failure")
  fails <- function(X, x0, ctx) stop(sentinel)
  fit$.model$distance$eval <- fit$.model$inference_distance$eval <- fails
  fit$.model$distance$eval_batch <- fit$.model$inference_distance$eval_batch <- NULL
  expect_identical(tryCatch(predict(fit, x0 = .25), error = identity), sentinel)
  expect_identical(tryCatch(vcov(fit), error = identity), sentinel)
})

test_that("valid current and candidate controls survive serialization and reuse", {
  for (method in c("exact", "asymp")) for (fit in .retained_control_fits(method)) {
    views <- list(fit)
    if (!is.null(fit$tuning_grid)) views <- c(views, list(candidate_fit(fit, 1L)))
    for (view in views) {
      restored <- unserialize(serialize(view, NULL))
      expect_identical(coef(restored), coef(view))
      if (!inherits(view, "rdistnn_gmm_result")) {
        expect_equal(as.numeric(predict(restored, X0 = restored$x0)), as.numeric(coef(restored)))
      }
      fields <- if (inherits(view, "rdistnn_tdnn_np_result")) c("s1", "s2", "dim") else c("k", "s")
      for (field in fields) {
        expect_equal(restored$meta[[field]], restored$.model[[field]])
        # Integer-valued doubles are also valid current control scalars.
        restored$meta[[field]] <- as.numeric(restored$meta[[field]])
        restored$.model[[field]] <- as.numeric(restored$.model[[field]])
      }
      expect_equal(coef(restored), coef(view))
    }
  }
})

test_that("reports require complete retained scientific state across estimator families", {
  fits <- .retained_report_fits()
  fields <- c(np = "y_vec", tdnn = "y_vec", cate = "psi", ratio = "signal",
    clate = "signal", plqr = "mD_mat")
  for (family in names(fields)) {
    fit <- fits[[family]]
    bad <- fit
    bad$.model[fields[[family]]] <- NULL
    routes <- c(.retained_report_routes(bad), list(print = function(x) capture.output(print(x))))
    for (route in routes) {
      error <- tryCatch(.with_localization_policy_warnings(route(bad)), error = identity)
      expect_s3_class(error, "rdistnn_state_error")
      if (inherits(error, "rdistnn_state_error")) {
        expect_identical(error$code, "state_invalid")
        expect_match(conditionMessage(error), "[Rr]efit")
      }
    }
    bad <- fit
    bad$.model[[fields[[family]]]] <- "malformed scientific state"
    for (route in routes) {
      error <- tryCatch(.with_localization_policy_warnings(route(bad)), error = identity)
      expect_s3_class(error, "rdistnn_state_error")
    }
  }
  for (field in c("signal_source", "created_with_version")) {
    bad <- fits$cate
    bad$meta[field] <- NULL
    for (route in .retained_report_routes(bad)) {
      expect_error(.with_localization_policy_warnings(route(bad)), class = "rdistnn_state_error")
    }
  }
})

test_that("retained model companions and controls reject malformed current values", {
  fits <- .retained_report_fits()
  cases <- list(
    np = list(y_vec = matrix(0, 8L), kind = "unknown", s = 9L, normalize = "TRUE"),
    tdnn = list(y_vec = rep(NA_real_, 8L), s1 = 4L, dim = 0L),
    cate = list(psi = matrix(0, 8L), k = 0L),
    ratio = list(row_ids = rep(1L, 8L), denominator_tol = -1, signal = matrix(0, 8L, 1L), s = 9L),
    clate = list(row_ids = NULL, denominator_tol = Inf, k = 0L),
    plqr = list(d_mat = matrix(0, 7L), mD_mat = matrix(Inf, 8L),
      theta_eval = "0", tau = 2, d_resid_mat = matrix(1, 8L),
      below_threshold = rep(2, 8L), density = rep(0, 8L), coef_names = character(),
      x_cols = 1L))
  for (family in names(cases)) for (field in names(cases[[family]])) {
    bad <- fits[[family]]
    bad$.model[field] <- list(cases[[family]][[field]])
    expect_error(coef(bad), class = "rdistnn_state_error")
  }
  for (field in c("d_resid_mat", "below_threshold", "density", "coef_names", "x_cols")) {
    bad <- fits$plqr
    bad$.model[field] <- NULL
    expect_error(summary(bad), class = "rdistnn_state_error")
  }
  for (field in c("row_ids", "denominator_tol")) {
    bad <- fits$ratio
    bad$.model[field] <- NULL
    expect_error(summary(bad), class = "rdistnn_state_error")
  }
  for (field in c("signal_source", "created_with_version")) {
    bad <- fits$cate
    bad$meta[[field]] <- c("invalid", "metadata")
    expect_error(coef(bad), class = "rdistnn_state_error")
  }
  for (field in c("theta_source", "p_d", "theta_se_method", "theta_meat_method", "moment_se_method")) {
    bad <- fits$plqr
    bad$meta[field] <- NULL
    expect_error(coef(bad), class = "rdistnn_state_error")
  }
  for (field in c("D", "one_sided")) {
    bad <- fits$clate
    bad$.model$iv[field] <- NULL
    expect_error(coef(bad), class = "rdistnn_state_error")
  }
  bad <- fits$clate
  bad$.model$iv$D[1L] <- 2
  expect_error(coef(bad), class = "rdistnn_state_error")
  bad <- fits$clate
  bad$.model$iv$one_sided <- "FALSE"
  expect_error(coef(bad), class = "rdistnn_state_error")
  bad <- fits$clate
  bad$.model$iv$propensity$raw <- NULL
  expect_error(coef(bad), class = "rdistnn_state_error")
  bad <- fits$clate
  bad$.model$iv$propensity$raw <- rep(Inf, 8L)
  expect_error(coef(bad), class = "rdistnn_state_error")
  disabled <- fits$plqr
  disabled["theta_se"] <- list(NULL)
  disabled$meta$theta_se_method <- "none"
  disabled$meta$theta_meat_method <- "none"
  expect_error(coef(disabled), class = "rdistnn_state_error")
})

test_that("current canonical report arrays reject malformed state on public consumers", {
  fits <- .retained_report_fits()
  for (family in names(fits)) {
    estimate <- if (family == "plqr") "theta" else if (family %in% c("ratio", "clate")) "theta" else if (family == "cate") "tau" else "mu"
    uncertainty <- if (family == "plqr") "theta_se" else "se"
    for (case in c("missing_estimate", "short_estimate", "nonfinite_estimate",
                   "missing_se", "short_se", "negative_se", "infinite_se",
                   "missing_targets", "short_targets", "nonfinite_targets", "target_count")) {
      bad <- fits[[family]]
      if (case == "missing_estimate") bad[estimate] <- NULL
      if (case == "short_estimate") bad[[estimate]] <- if (family == "plqr") matrix(0, 1L, 1L) else 0
      if (case == "nonfinite_estimate") bad[[estimate]][1L] <- Inf
      if (case == "missing_se") bad[uncertainty] <- NULL
      if (case == "short_se") bad[[uncertainty]] <- if (family == "plqr") matrix(1, 1L, 1L) else 1
      if (case == "negative_se") bad[[uncertainty]][1L] <- -2
      if (case == "infinite_se") bad[[uncertainty]][1L] <- Inf
      if (case == "missing_targets") bad["x0"] <- NULL
      if (case == "short_targets") bad$x0 <- bad$x0[1L, , drop = FALSE]
      if (case == "nonfinite_targets") bad$x0[1L] <- NA_real_
      if (case == "target_count") bad$meta$m <- 1L
      before <- serialize(bad, NULL)
      for (route in .retained_report_routes(bad)) {
        error <- tryCatch({ .with_localization_policy_warnings(route(bad)); NULL }, error = identity)
        expect_true(inherits(error, "error") && grepl("[Rr]efit", conditionMessage(error)),
          info = paste(family, case))
      }
      .expect_identical_fast(serialize(bad, NULL), before, info = paste(family, case))
    }
  }
})

test_that("valid current arrays survive serialization, repeated targets and disabled uncertainty", {
  .with_partial_result_warning({
  for (se in c(TRUE, FALSE)) for (fit in .retained_report_fits(se)) {
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(coef(restored), coef(fit))
    expect_equal(summary(restored), summary(fit))
    expect_equal(.with_localization_policy_warnings(plot_data_estimates(restored)), .with_localization_policy_warnings(plot_data_estimates(fit)))
    expect_equal(.with_localization_policy_warnings(predict(restored, x0 = .25)), .with_localization_policy_warnings(predict(fit, x0 = .25)))
    if (!is.null(fit$tuning_grid)) {
      candidate <- .with_localization_policy_warnings(candidate_fit(restored, 1L))
      reference <- .with_localization_policy_warnings(candidate_fit(fit, 1L))
      expect_equal(coef(candidate), coef(reference))
      expect_equal(summary(candidate), summary(reference))
      expect_identical(candidate$meta$candidate_id, reference$meta$candidate_id)
      numerical <- names(reference$.model)[vapply(reference$.model, is.numeric, logical(1L))]
      .expect_identical_fast(candidate$.model[numerical], reference$.model[numerical])
      expect_s3_class(summary(candidate), class(summary(fit)))
      expect_equal(.with_localization_policy_warnings(predict(candidate, x0 = .25)),
        .with_localization_policy_warnings(predict(candidate_fit(fit, 1L), x0 = .25)))
    }
  }
  fit <- affine_ratio_dnn(X = matrix(1:4), X0 = matrix(c(0, 5)),
    numerator = rep(1, 4), denominator = c(0, 0, 1, 1), s = 4, se = FALSE)
  expect_true(is.na(coef(fit)[1L]))
  expect_s3_class(summary(fit), "rdistnn_summary_affine_ratio_result")
  expect_equal(.with_localization_policy_warnings(predict(fit, X0 = fit$x0)), unname(coef(fit)))
  failed <- affine_ratio_dnn(X = matrix(1:4), x0 = 0,
    numerator = rep(1, 4), denominator = rep(0, 4), s = 2)
  expect_true(all(is.na(summary(failed)$coefficients$std.error)))
  })
})

test_that("paired and PLQR covariance arrays retain target and component dimensions", {
  fits <- .retained_report_fits()[c("ratio", "clate", "plqr")]
  for (family in names(fits)) {
    fields <- if (family == "plqr") "theta_vcov" else c("theta_vcov", "component_vcov")
    for (field in fields) for (case in c("short_list", "wrong_matrix", "negative_diagonal")) {
      bad <- fits[[family]]
      if (case == "short_list") bad[[field]] <- bad[[field]][1L]
      if (case == "wrong_matrix") bad[[field]][[1L]] <- matrix(0, 3L, 3L)
      if (case == "negative_diagonal") bad[[field]][[1L]][1L, 1L] <- -1
      for (route in .retained_report_routes(bad)) {
        error <- tryCatch({ .with_localization_policy_warnings(route(bad)); NULL }, error = identity)
        expect_true(inherits(error, "error") && grepl("[Rr]efit", conditionMessage(error)),
          info = paste(family, field, case))
      }
    }
  }
})


test_that("moment-only PLQR candidate views keep explicit disabled uncertainty", {
  tuning <- new_tuning("retained_moment_candidate", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L)))
  fit <- plqr_dnn(X = matrix(1:8), x0 = 0, Y = rep(1, 8), D = matrix(1:8),
    mD = matrix(0, 8, 1), g = rep(0, 8), theta_eval = 0, se = FALSE,
    tuning = tuning, tuning_mode = "grid", return = character())
  candidate <- candidate_fit(fit, 1L)
  expect_true("moment_se" %in% names(candidate))
  expect_null(candidate$moment_se)
  expect_s3_class(summary(candidate), "rdistnn_summary_plqr_result")
})

test_that("scalar ratio candidate plot tables reject duplicate target rows", {
  fits <- .retained_report_fits()[c("ratio", "clate")]
  for (fit in fits) {
    bad <- fit
    bad$tuning_grid$estimates <- rbind(bad$tuning_grid$estimates,
      bad$tuning_grid$estimates[1L, , drop = FALSE])
    expect_error(.with_localization_policy_warnings(plot_data_estimates(bad, candidates = "all")), "duplicate target_id")
    expect_error(.with_localization_policy_warnings(plot_data_estimates(bad, candidates = bad$tuning_grid$candidates$candidate_id[[1L]])), "duplicate target_id")
    expect_equal(nrow(.with_localization_policy_warnings(plot_data_estimates(fit, candidates = "all"))),
      nrow(fit$tuning_grid$candidates) * nrow(fit$x0))
  }
})

test_that("paired uncertainty arrays agree with the retained covariance mode", {
  .with_partial_result_warning({
  fields <- c("se", "component_vcov", "theta_vcov", "theta_omega")
  for (fit in .retained_report_fits()[c("ratio", "clate")]) {
    expect_identical(fit$meta$covariance, "paired_frozen_score")
    expect_true(any(fit$status$wald_available))
    corruptions <- lapply(fields, function(field) {
      bad <- fit
      bad[field] <- list(NULL)
      bad
    })
    jointly_null <- fit
    jointly_null[fields] <- rep(list(NULL), length(fields))
    jointly_null$status$wald_available <- FALSE
    corruptions <- c(corruptions, list(jointly_null))
    for (bad in corruptions) for (route in .retained_report_routes(bad)) {
      expect_error(.with_localization_policy_warnings(route(bad)), "[Rr]efit")
    }
    for (mode in list(NULL, "", "unknown", c("none", "paired_frozen_score"))) {
      bad <- fit
      bad$meta["covariance"] <- list(mode)
      expect_error(summary(bad), "meta\\$covariance.*[Rr]efit")
    }
    disabled <- fit
    disabled$meta$covariance <- "none"
    expect_error(summary(disabled), "[Rr]efit")
  }
  for (fit in .retained_report_fits(FALSE)[c("ratio", "clate")]) {
    expect_identical(fit$meta$covariance, "none")
    expect_true(all(vapply(fit[fields], is.null, logical(1L))))
    expect_s3_class(summary(fit), "rdistnn_summary_affine_ratio_result")
    bad <- fit
    bad$status$wald_available[[1L]] <- TRUE
    expect_error(summary(bad), "[Rr]efit")
  }
  unavailable <- affine_ratio_dnn(X = matrix(1:4), x0 = 0,
    numerator = c(1, 2, 3, 4) * 1e-200, denominator = rep(1, 4), s = 2)
  expect_identical(unavailable$status$covariance_status, "numerical_failure")
  expect_true(is.na(unavailable$se))
  expect_identical(unavailable$meta$covariance, "paired_frozen_score")
  expect_s3_class(summary(unavailable), "rdistnn_summary_affine_ratio_result")
  expect_true(all(is.na(confint(unavailable))))
  })
})

test_that("ordinary scalar uncertainty agrees with retained jackknife mode", {
  for (fit in .retained_report_fits()[c("np", "cate", "tdnn")]) {
    expect_identical(fit$meta$jackknife, "fast")
    bad <- fit
    bad["se"] <- list(NULL)
    for (route in .retained_report_routes(bad)) {
      expect_error(.with_localization_policy_warnings(route(bad)), "[Rr]efit")
    }
    disabled <- fit
    disabled$meta$jackknife <- "none"
    expect_error(summary(disabled), "[Rr]efit")
    for (mode in list(NULL, "", "unknown", c("fast", "none"))) {
      bad <- fit
      bad$meta["jackknife"] <- list(mode)
      expect_error(summary(bad), "meta\\$jackknife.*[Rr]efit")
    }
  }
  for (fit in .retained_report_fits(FALSE)[c("np", "cate", "tdnn")]) {
    expect_identical(fit$meta$jackknife, "none")
    expect_null(fit$se)
    expect_s3_class(summary(fit), if (inherits(fit, "rdistnn_cate_result")) "rdistnn_summary_cate_result" else "rdistnn_summary_np_result")
  }
})

test_that("paired reports reject missing and inconsistent retained jackknife modes", {
  for (se in c(TRUE, FALSE)) for (fit in .retained_report_fits(se)[c("ratio", "clate")]) {
    modes <- c(list(NULL, "", "unknown", c("fast", "none")),
      list(if (se) "none" else "fast"))
    for (mode in modes) {
      bad <- fit
      bad$meta["jackknife"] <- list(mode)
      for (route in .retained_report_routes(bad)) {
        expect_error(.with_localization_policy_warnings(route(bad)),
          class = "rdistnn_state_error")
      }
    }
  }
})

test_that("paired retained weight methods reject missing, invalid and inconsistent controls", {
  fits <- .retained_report_fits(method = "asymp")[c("ratio", "clate")]
  exact <- .retained_report_fits(method = "exact")[c("ratio", "clate")]
  for (family in names(fits)) {
    fit <- fits[[family]]
    # This fixture distinguishes the former NULL-to-exact prediction fallback.
    expect_gt(max(abs(predict(fit, X0 = fit$x0) - predict(exact[[family]], X0 = fit$x0))), .01)
    cases <- c("missing", "null", "na", "vector", "invalid", "mismatch")
    for (location in c(".model", "meta")) for (case in cases) {
      bad <- fit
      if (case == "missing") bad[[location]]["method"] <- NULL
      if (case == "null") bad[[location]]["method"] <- list(NULL)
      if (case == "na") bad[[location]]$method <- NA_character_
      if (case == "vector") bad[[location]]$method <- c("exact", "asymp")
      if (case == "invalid") bad[[location]]$method <- "unrecognized"
      if (case == "mismatch") bad[[location]]$method <- "exact"
      before <- serialize(bad, NULL)
      routes <- c(.retained_report_routes(bad), list(
        prediction_se = function(x) predict(x, X0 = x$x0, se.fit = TRUE),
        prediction_interval = function(x) predict(x, X0 = x$x0, interval = "confidence")))
      for (route in routes) {
        error <- tryCatch({ .with_localization_policy_warnings(route(bad)); NULL }, error = identity)
        expect_true(inherits(error, "error") && grepl("[Rr]efit", conditionMessage(error)),
          info = paste(family, location, case))
      }
      .expect_identical_fast(serialize(bad, NULL), before, info = paste(family, location, case))
    }
  }
})

test_that("paired exact and asymptotic controls survive serialized point and SE prediction", {
  for (method in c("exact", "asymp")) for (se in c(TRUE, FALSE)) {
    fits <- .retained_report_fits(se = se, method = method)[c("ratio", "clate")]
    for (fit in fits) {
      expect_identical(fit$.model$method, method)
      expect_identical(fit$meta$method, method)
      restored <- unserialize(serialize(fit, NULL))
      point <- predict(fit, X0 = fit$x0)
      uncertainty <- predict(fit, X0 = fit$x0, se.fit = TRUE)
      expect_equal(point, unname(coef(fit)))
      expect_identical(predict(restored, X0 = fit$x0), point)
      expect_identical(predict(restored, X0 = fit$x0, se.fit = TRUE), uncertainty)
      if (se) expect_equal(uncertainty$se.fit, fit$se)
      candidate <- candidate_fit(restored, 1L)
      expect_identical(candidate$.model$method, method)
      expect_identical(candidate$meta$method, method)
      expect_equal(predict(candidate, X0 = fit$x0, se.fit = TRUE),
        predict(candidate_fit(fit, 1L), X0 = fit$x0, se.fit = TRUE))
    }
  }
})
