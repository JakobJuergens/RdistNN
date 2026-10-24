.candidate_retained_state_fixtures <- function(extras = character()) {
  X <- matrix(c(0, 0, 0, 1, 2, 2, 3, 4))
  X0 <- matrix(c(0, 2))
  y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  tuning <- new_tuning("candidate_state_validation",
    eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 3L, 5L)))
  common <- list(X = X, X0 = X0, return = extras)
  constructors <- list(
    np = function(...) np_dnn(Y = y, ...),
    cate = function(...) cate_dnn(pseudo_outcome = y, ...),
    plqr = function(...) plqr_dnn(Y = y, D = matrix(seq_len(8L)),
      mD = matrix(0, 8L, 1L), g = rep(0, 8L), theta_eval = .3, ...),
    ratio = function(...) affine_ratio_dnn(numerator = y,
      denominator = c(1, 3, .5, 2, 4, 2, 1, 3), ...)
  )
  fits <- lapply(constructors, function(fn) {
    .with_distance_tie_warning(do.call(fn,
      c(common, list(tuning = tuning, tuning_mode = "grid"))))
  })
  list(fits = fits, constructors = constructors, common = common)
}

test_that("candidate-state fixtures preserve unrelated estimator warnings", {
  for (estimator in c("np_dnn", "cate_dnn", "plqr_dnn", "affine_ratio_dnn")) {
    real_estimator <- get(estimator, mode = "function")
    sentinel <- simpleWarning("unrelated candidate-state fixture warning")
    fixture <- .candidate_retained_state_fixtures
    environment(fixture) <- list2env(stats::setNames(list(function(...) {
      warning(sentinel)
      real_estimator(...)
    }), estimator), parent = environment(fixture))
    seen <- list()
    withCallingHandlers(fixture(), warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    expect_equal(length(seen), 1L, info = estimator)
    if (length(seen)) expect_identical(seen[[1L]], sentinel, info = estimator)
  }
})

test_that("retained candidate scalar corruption is a state error with refit guidance", {
  fit <- .candidate_retained_state_fixtures()$fits$np
  for (field in c("k", "s", "normalize")) {
    broken <- fit
    broken$.model[[field]] <- if (field == "normalize") NA else 0L
    err <- tryCatch(.with_localization_policy_warnings(candidate_fit(broken, 1L)), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    expect_match(conditionMessage(err), "Refit", fixed = TRUE)
  }
})

test_that("candidate payloads require nonempty unique candidate IDs", {
  fits <- .candidate_retained_state_fixtures()$fits
  for (fit in fits) {
    for (invalid in list(NULL, rep(NA_character_, 4L), rep("", 4L), rep("same", 4L))) {
      bad <- fit
      bad$tuning_grid$candidates$candidate_id <- invalid
      err <- tryCatch(.with_localization_policy_warnings(candidate_fit(bad, 1L)), error = identity)
      expect_s3_class(err, "rdistnn_state_error")
      expect_identical(err$code, "state_invalid")
      expect_match(conditionMessage(err), "candidate_id", fixed = TRUE)
      expect_match(conditionMessage(err), "Refit", fixed = TRUE)
    }
  }
})

test_that("candidate consumers reject incomplete retained selection schemas", {
  fit <- .candidate_retained_state_fixtures()$fits$np
  consumers <- list(function(x) candidate_fit(x, 1L),
    function(x) diagnostics(x, candidates = "all"),
    function(x) plot_data_estimates(x, candidates = "all"))
  for (field in c("status", "reason", "is_default", "feasible", "computed_key", "s", "lambda")) {
    bad <- fit
    bad$tuning_grid$candidates[[field]] <- NULL
    for (consume in consumers) {
      err <- tryCatch(.with_localization_policy_warnings(consume(bad)), error = identity)
      expect_s3_class(err, "rdistnn_state_error")
      if (inherits(err, "condition")) {
        expect_match(conditionMessage(err), field, fixed = TRUE)
        expect_match(conditionMessage(err), "Refit", fixed = TRUE)
      }
    }
  }
})

test_that("candidate consumers reject malformed stored scientific values before coercion", {
  fits <- .candidate_retained_state_fixtures()$fits[c("np", "cate", "plqr")]
  consumers <- list(function(x) candidate_fit(x, 1L),
    function(x) plot_data_estimates(x, candidates = "all"))
  for (fit in fits) {
    for (field in c("estimate", "se")) {
      invalid <- if (field == "estimate") list("not numeric", Inf, NaN, NA_real_) else
        list("not numeric", Inf, NaN, NA_real_, -1)
      for (value in invalid) {
        bad <- fit
        bad$tuning_grid$estimates[[field]] <- rep(value, nrow(bad$tuning_grid$estimates))
        for (consume in consumers) {
          err <- tryCatch(.with_localization_policy_warnings(consume(bad)), error = identity)
          expect_s3_class(err, "rdistnn_state_error")
          if (inherits(err, "condition")) {
            expect_match(conditionMessage(err), field, fixed = TRUE)
            expect_match(conditionMessage(err), "Refit", fixed = TRUE)
          }
        }
      }
    }
    bad <- fit
    bad$tuning_grid$estimates$component <- "unknown component"
    err <- tryCatch(.with_localization_policy_warnings(candidate_fit(bad, 1L)), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    if (inherits(err, "condition")) expect_match(conditionMessage(err), "component", fixed = TRUE)
  }
  point <- fits$np
  point <- .with_localization_policy_warnings(np_dnn(X = point$.model$x_mat,
    X0 = point$x0, Y = point$.model$y_vec, se = FALSE,
    tuning = new_tuning("point_grid", eval = function(...) list(s = 3L),
      eval_grid = function(...) data.frame(s = c(2L, 3L))), tuning_mode = "grid"))
  point$tuning_grid$estimates$se <- 1
  err <- tryCatch(.with_localization_policy_warnings(candidate_fit(point, 1L)), error = identity)
  expect_s3_class(err, "rdistnn_state_error")
  if (inherits(err, "condition")) expect_match(conditionMessage(err), "se", fixed = TRUE)
})

test_that("retained long estimates validate every consumed column before selection", {
  fits <- .candidate_retained_state_fixtures()$fits[c("np", "cate", "ratio")]
  consumers <- list(coef, function(x) candidate_fit(x, 1L),
    function(x) plot_data_estimates(x, candidates = "all"))
  invalid <- list(lambda = list("bad", TRUE, 1i),
    status = list(TRUE, NA_character_, "failed"), reason = list(1L, NA_character_, "different reason"))
  for (fit in fits) for (field in names(invalid)) {
    for (value in c(list(NULL), invalid[[field]])) {
      bad <- fit
      bad$tuning_grid$estimates[[field]] <- if (is.null(value)) NULL else
        rep(value, nrow(bad$tuning_grid$estimates))
      for (consume in consumers) {
        error <- tryCatch(.with_localization_policy_warnings(consume(bad)), error = identity)
        expect_s3_class(error, "rdistnn_state_error")
        if (inherits(error, "condition")) {
          expect_match(conditionMessage(error), field, fixed = TRUE)
          expect_match(conditionMessage(error), "[Rr]efit")
        }
      }
    }
  }
})

test_that("retained candidate diagnostics validate types and availability before selection", {
  fit <- .candidate_retained_state_fixtures()$fits$np
  rounded <- fit
  rounded$tuning_grid$diagnostics$top5_weight_share <- 1 + 4 * .Machine$double.eps
  expect_no_error(.with_localization_policy_warnings(diagnostics(rounded, candidates = "all")))
  invalid <- list(lambda = list("bad", TRUE, 1i), degenerate = list("yes", 1),
    status = list(TRUE, NA_character_, "failed"), reason = list(1L, NA_character_, "different reason"),
    ess = list("bad", -1), max_weight = list("bad", -1), top5_weight_share = list("bad", -1, 2))
  for (field in names(invalid)) for (value in invalid[[field]]) {
    bad <- fit
    bad$tuning_grid$diagnostics[[field]] <- rep(value, nrow(bad$tuning_grid$diagnostics))
    for (consume in list(coef, function(x) diagnostics(x, candidates = "all"))) {
      error <- tryCatch(.with_localization_policy_warnings(consume(bad)), error = identity)
      expect_s3_class(error, "rdistnn_state_error")
      if (inherits(error, "condition")) {
        expect_match(conditionMessage(error), field, fixed = TRUE)
        expect_match(conditionMessage(error), "[Rr]efit")
      }
    }
  }
})

test_that("current numeric candidate labels and unavailable diagnostic rows survive serialization", {
  fixture <- .candidate_retained_state_fixtures()
  labels <- c(NA_real_, NaN, Inf, -Inf, 1)
  tuning <- new_tuning("numeric_candidate_labels", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(lambda = labels, s = c(2L, 3L, 4L, 5L, 9L)))
  for (fn in fixture$constructors[c("np", "cate", "ratio")]) {
    fit <- .with_test_warning_families(.with_localization_policy_warnings(do.call(fn,
      c(fixture$common, list(tuning = tuning, tuning_mode = "grid")))), "infeasible_candidates")
    saved <- unserialize(serialize(fit, NULL))
    expect_identical(saved$tuning_grid$estimates$lambda, fit$tuning_grid$estimates$lambda)
    expect_true(any(is.na(saved$tuning_grid$diagnostics$degenerate)))
    expect_equal(.with_localization_policy_warnings(coef(saved)),
      .with_localization_policy_warnings(coef(fit)))
    expect_s3_class(.with_localization_policy_warnings(diagnostics(saved, candidates = "all")),
      "rdistnn_diagnostics")
    expect_s3_class(.with_localization_policy_warnings(plot_data_estimates(saved, candidates = "all")),
      "data.frame")
    expect_no_error(.with_localization_policy_warnings(candidate_fit(saved, 1L)))
  }
})

test_that("retained grid defaults must succeed while fresh sensitivity anchors can be infeasible", {
  fit <- .candidate_retained_state_fixtures()$fits$np
  bad <- fit
  id <- bad$tuning_grid$candidates$candidate_id[bad$tuning_grid$candidates$is_default]
  for (field in c("candidates", "estimates", "diagnostics")) {
    selected <- bad$tuning_grid[[field]]$candidate_id == id
    bad$tuning_grid[[field]]$status[selected] <- "failed"
    bad$tuning_grid[[field]]$reason[selected] <- "failure"
  }
  for (consume in list(coef, function(x) diagnostics(x, candidates = "all"))) {
    error <- expect_error(.with_localization_policy_warnings(consume(bad)), class = "rdistnn_state_error")
    expect_identical(error$code, "state_invalid")
    expect_match(conditionMessage(error), "[Rr]efit")
  }
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    scalar <- np_dnn(X = matrix(1:8), x0 = .25, Y = (1:8)^2, s = 8L, se = FALSE)
    expect_warning(panel <- plot_tuning(scalar, se = TRUE), class = "rdistnn_infeasible_candidates_warning")
    expect_s3_class(panel, "ggplot")
  }
})

test_that("candidate selection errors bound long lists and escape literal names", {
  fit <- .candidate_retained_state_fixtures()$fits$np
  unknown <- paste0("unknown_", seq_len(200L))
  unknown[1L] <- paste0("line\nbreak_", strrep("x", 150L))
  for (consume in list(function() candidate_fits(fit, unknown),
      function() plot_data_estimates(fit, candidates = unknown))) {
    err <- tryCatch(.with_localization_policy_warnings(consume()), error = identity)
    expect_s3_class(err, "rdistnn_input_error")
    expect_identical(err$code, "argument_invalid")
    expect_lt(nchar(conditionMessage(err)), 800L)
    expect_match(conditionMessage(err), "195 more", fixed = TRUE)
    expect_false(grepl("\n", conditionMessage(err), fixed = TRUE))
  }
})

test_that("unchanged candidate schemas are validated once per public reuse call", {
  fit <- .candidate_retained_state_fixtures()$fits$np
  checks <- 0L
  original <- RdistNN:::.dnn_grid_check_candidates
  local_mocked_bindings(.dnn_grid_check_candidates = function(candidates, retained = FALSE) {
    if (retained) checks <<- checks + 1L
    original(candidates, retained)
  }, .package = "RdistNN")
  for (i in 1:2) {
    checks <- 0L
    views <- .with_localization_policy_warnings(candidate_fits(fit, 1:3))
    expect_length(views, 3L)
    expect_identical(checks, 1L)
  }
  changed <- fit
  changed$tuning_grid$candidates$status <- NULL
  expect_error(.with_localization_policy_warnings(candidate_fit(changed, 1L)),
    class = "rdistnn_state_error")
})

test_that("the retained default candidate matches the fitted order", {
  fit <- .candidate_retained_state_fixtures()$fits$np
  fit$tuning_grid$candidates$is_default <- c(TRUE, FALSE, FALSE, FALSE)
  for (consume in list(function(x) candidate_fit(x, 2L),
      function(x) diagnostics(x, candidates = "all"),
      function(x) plot_data_estimates(x, candidates = "all"))) {
    err <- tryCatch(.with_localization_policy_warnings(consume(fit)), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    if (inherits(err, "condition")) expect_match(conditionMessage(err), "default", fixed = TRUE)
  }
})

test_that("paired candidate status facts are validated before interval indexing", {
  fit <- .candidate_retained_state_fixtures()$fits$ratio
  expect_true(any(fit$tuning_grid$estimates$wald_available))
  corruptions <- list(
    function(x) { x$wald_available <- "yes"; x },
    function(x) { x$wald_available <- NA; x },
    function(x) { x$wald_available <- 1; x },
    function(x) { x$wald_available <- FALSE; x },
    function(x) { x$ratio_status <- NA_character_; x },
    function(x) { x$ratio_status <- "zero_denominator"; x },
    function(x) { x$covariance_status <- "not_computed"; x }
  )
  for (change in corruptions) {
    bad <- fit
    bad$tuning_grid$estimates <- change(bad$tuning_grid$estimates)
    err <- tryCatch(.with_localization_policy_warnings(plot_data_estimates(bad, candidates = "all")),
      error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    if (inherits(err, "condition")) expect_match(conditionMessage(err), "Refit", fixed = TRUE)
  }
})

test_that("candidate estimate targets are validated before integer conversion", {
  # Ratios recompute paired estimates from their frozen signals; these three
  # adapters consume the stored estimate rows.
  fits <- .candidate_retained_state_fixtures()$fits[c("np", "cate", "plqr")]
  for (fit in fits) {
    for (change in list(function(x) x + .25, function(x) as.character(x),
                       function(x) rep(Inf, length(x)), function(x) as.complex(x) + 1i)) {
      bad <- fit
      bad$tuning_grid$estimates$target_id <- change(bad$tuning_grid$estimates$target_id)
      error <- expect_error(suppressWarnings(candidate_fit(bad, 1L)), class = "rdistnn_state_error")
      expect_identical(error$code, "state_invalid")
      expect_match(conditionMessage(error), "target_id", fixed = TRUE)
    }
  }
})

test_that("selected successful candidate orders are validated before conversion", {
  fits <- .candidate_retained_state_fixtures()$fits
  for (fit in fits) {
    for (invalid in list(3.5, Inf, NaN, NA_real_, 2^31, -2^31, 0, -1, "3", 3 + 1i)) {
      bad <- fit
      bad$tuning_grid$candidates$s[1L] <- invalid
      expect_error(suppressWarnings(candidate_fit(bad, 1L)),
        "candidates.*s.*positive integer")
    }
    # Unavailable rows are retained for grid reporting; extracting a separate
    # successful row must not try to consume their unavailable order.
    unavailable <- fit
    unavailable$tuning_grid$candidates$s[4L] <- NA_integer_
    unavailable$tuning_grid$candidates$status[4L] <- "infeasible"
    unavailable$tuning_grid$candidates$feasible[4L] <- FALSE
    unavailable$tuning_grid$candidates$reason[4L] <- "unavailable order"
    id <- unavailable$tuning_grid$candidates$candidate_id[4L]
    for (field in c("estimates", "diagnostics")) {
      selected <- unavailable$tuning_grid[[field]]$candidate_id == id
      unavailable$tuning_grid[[field]]$status[selected] <- "infeasible"
      unavailable$tuning_grid[[field]]$reason[selected] <- "unavailable order"
    }
    expect_equal(suppressWarnings(candidate_fit(unavailable, 1L)),
      suppressWarnings(candidate_fit(fit, 1L)))
    expect_error(suppressWarnings(candidate_fit(unavailable, 4L)), "infeasible.*unavailable order")
    whole <- fit
    whole$tuning_grid$candidates$s <- as.numeric(whole$tuning_grid$candidates$s)
    expect_identical(suppressWarnings(candidate_fit(whole, 1L))$meta$s, 2L)
  }
})

test_that("candidate dimension metadata rejects integer overflow cleanly", {
  fits <- .candidate_retained_state_fixtures()$fits
  for (fit in fits) for (field in c("n", "m")) {
    for (invalid in c(2^31, -2^31, 1e20, -1e20)) {
      bad <- fit
      bad$meta[[field]] <- invalid
      expect_error(suppressWarnings(candidate_fit(bad, 1L)),
        paste0("object\\$meta\\$", field, ".*integer scalar"))
    }
  }
})

test_that("thin candidate views validate cached rank shapes and permutations", {
  fits <- .candidate_retained_state_fixtures()$fits
  for (fit in fits) {
    for (invalid in list(matrix(rep(1:7, 2L), 7L, 2L), matrix(1:8),
                         matrix(rep(c(1L, 1L, 3:8), 2L), 8L, 2L))) {
      bad <- fit
      bad$.model$ranks <- invalid
      expect_error(suppressWarnings(candidate_fit(bad, 1L)), "[Rr]anks|permutation")
    }
  }
})

test_that("serialized thin candidate views preserve valid duplicate-order rows", {
  fixture <- .candidate_retained_state_fixtures()
  value <- function(fit) {
    if (inherits(fit, "rdistnn_plqr_result")) fit$moment else stats::coef(fit)
  }
  for (kind in names(fixture$fits)) {
    fit <- unserialize(serialize(fixture$fits[[kind]], NULL))
    candidates <- fit$tuning_grid$candidates
    expect_identical(anyDuplicated(candidates$candidate_id), 0L)
    expect_gt(anyDuplicated(candidates$computed_key), 0L)
    for (i in seq_len(nrow(candidates))) {
      view <- suppressWarnings(candidate_fit(fit, i))
      scalar <- suppressWarnings(do.call(fixture$constructors[[kind]],
        c(fixture$common, list(s = candidates$s[i]))))
      expect_identical(class(view), class(fit))
      expect_identical(view$meta$candidate_id, candidates$candidate_id[i])
      expect_identical(view$meta$s, candidates$s[i])
      expect_null(view$W)
      expect_equal(suppressWarnings(value(view)), suppressWarnings(value(scalar)), tolerance = 1e-12)
      expect_equal(suppressWarnings(predict(view, X0 = fixture$common$X0)),
        suppressWarnings(predict(scalar, X0 = fixture$common$X0)), tolerance = 1e-12)
    }
  }
})

test_that("thin candidate views require valid retained controls independently of extras", {
  for (extras in list(character(), "W", c("ess", "neighbors"))) {
    fixture <- .candidate_retained_state_fixtures(extras)
    for (fit in fixture$fits[c("np", "cate", "plqr")]) {
      for (field in c("x_mat", "method", "k", "s", "normalize")) {
        bad <- fit
        bad$.model[[field]] <- NULL
        expect_error(suppressWarnings(candidate_fit(bad, 1L)), "Refit|refit")
      }
      for (field in c("method", "k", "s")) {
        bad <- fit
        bad$meta[[field]] <- switch(field, method = "asymp", k = 2L, s = 5L)
        expect_error(suppressWarnings(candidate_fit(bad, 1L)), "inconsistent|match")
      }
    }
  }
})

test_that("candidate order identities cannot reuse stored rows for another order", {
  fits <- .candidate_retained_state_fixtures()$fits
  for (fit in fits) {
    bad <- fit
    bad$tuning_grid$candidates$s[1L] <- 7L
    expect_error(suppressWarnings(candidate_fit(bad, 1L)), "computed_key|inconsistent")
    bad <- fit
    bad$tuning_grid$candidates$s[1L] <- 9L
    expect_error(suppressWarnings(candidate_fit(bad, 1L)), "s.*n|inconsistent")
    bad <- fit
    id <- bad$tuning_grid$candidates$candidate_id[1L]
    bad$tuning_grid$estimates$s[bad$tuning_grid$estimates$candidate_id == id] <- 7L
    expect_error(suppressWarnings(candidate_fit(bad, 1L)), "order|inconsistent|match")
  }
})

test_that("candidate display consumers validate the same retained identities as extraction", {
  fits <- .candidate_retained_state_fixtures("W")$fits
  consumers <- list(
    function(fit) diagnostics(fit, candidates = "all"),
    function(fit) plot_data_estimates(fit, candidates = "all")
  )
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    consumers <- c(consumers, list(
      function(fit) plot_estimates(fit, candidates = "all"),
      function(fit) plot_weights(fit, candidates = "all")))
  }
  for (fit in fits) {
    corruptions <- list(
      function(x) { x$tuning_grid$candidates$s[1L] <- 7L; x },
      function(x) { x$tuning_grid$estimates$target_id[1L] <- 1.5; x },
      function(x) { x$meta$k <- 2L; x },
      function(x) { x$tuning_grid$candidates$candidate_id[2L] <- x$tuning_grid$candidates$candidate_id[1L]; x }
    )
    for (change in corruptions) {
      bad <- change(fit)
      for (consume in consumers) {
        expect_error(suppressWarnings(consume(bad)), "inconsistent|integer.*target_id|unique|match")
      }
    }
    unavailable <- fit
    unavailable$tuning_grid$candidates$s[4L] <- NA_integer_
    unavailable$tuning_grid$candidates$status[4L] <- "infeasible"
    unavailable$tuning_grid$candidates$feasible[4L] <- FALSE
    unavailable$tuning_grid$candidates$reason[4L] <- "unavailable order"
    id <- unavailable$tuning_grid$candidates$candidate_id[4L]
    for (field in c("estimates", "diagnostics")) {
      selected <- unavailable$tuning_grid[[field]]$candidate_id == id
      unavailable$tuning_grid[[field]]$status[selected] <- "infeasible"
      unavailable$tuning_grid[[field]]$reason[selected] <- "unavailable order"
    }
    for (consume in consumers) expect_no_error(suppressWarnings(consume(unavailable)))
  }
})

test_that("candidate diagnostics validate their recorded computation identity", {
  fits <- .candidate_retained_state_fixtures()$fits
  for (fit in fits) {
    for (field in c("s", "computed_key")) {
      bad <- fit
      bad$tuning_grid$diagnostics[[field]][1L] <- if (field == "s") 7L else "s_7"
      expect_error(suppressWarnings(diagnostics(bad, candidates = "all")),
        "diagnostics.*inconsistent.*[Rr]efit")
    }
  }
})

test_that("candidate views require current tuning metadata instead of repairing it", {
  fits <- .candidate_retained_state_fixtures(c("W", "tuning"))$fits
  changes <- list(
    function(x) { x$meta$tuning <- NULL; x },
    function(x) { x$meta$tuning <- NULL; x$tuning <- NULL; x },
    function(x) { names(x$meta)[names(x$meta) == "tuning"] <- "tuning_old"; x },
    function(x) { x$meta$tuning$default <- NULL; x },
    function(x) { x$meta$tuning$default$meta <- NULL; x },
    function(x) { x$meta$tuning$default$meta$name <- NULL; x },
    function(x) { x$meta$tuning$default$meta$n <- NULL; x },
    function(x) { x$meta$tuning$default$s <- 4L; x },
    function(x) { x$meta$tuning$default$feasibility <- NULL; x },
    function(x) { x$meta$tuning$default$feasibility$k_le_s <- FALSE; x },
    function(x) { x$meta$tuning$source <- NULL; x },
    function(x) { x$meta$tuning$grid <- NULL; x },
    function(x) { x$meta$tuning$grid$meta <- NULL; x }
  )
  for (fit in fits) for (change in changes) {
    bad <- change(fit)
    err <- tryCatch(.with_localization_policy_warnings(candidate_fit(bad, 1L)), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    if (inherits(err, "condition")) expect_match(conditionMessage(err), "Refit", fixed = TRUE)
  }
})

test_that("current tuning and optional manual tuning survive supported reuse", {
  fits <- .candidate_retained_state_fixtures(c("W", "tuning"))$fits
  for (fit in fits) {
    saved <- unserialize(serialize(fit, NULL))
    view <- .with_localization_policy_warnings(candidate_fit(saved, 1L))
    expect_identical(view$meta$tuning$default$s, 2L)
    expect_identical(view$meta$tuning$default$meta$name, fit$meta$tuning$default$meta$name)
    expect_identical(view$meta$tuning$default$meta$candidate_id,
      fit$tuning_grid$candidates$candidate_id[1L])
    expect_identical(view$tuning, view$meta$tuning)
  }
  manual <- np_dnn(X = matrix(1:8), x0 = .25, Y = (1:8)^2, s = 3L,
    se = FALSE, return = c("W", "tuning"))
  expect_null(manual$meta[["tuning", exact = TRUE]])
  expect_true(is.list(manual[["tuning", exact = TRUE]]))
  expect_identical(diagnostics(manual)$tuning_source, "manual")
})
