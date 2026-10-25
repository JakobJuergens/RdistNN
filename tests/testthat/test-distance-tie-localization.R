# Independent complete-subset reference: uniform choice at the kth represented
# distance boundary, with no production rank tables or group helpers.
.distance_subset_weights <- function(dx, k, s) {
  dx <- as.matrix(dx)
  subsets <- utils::combn(nrow(dx), s, simplify = FALSE)
  out <- matrix(0, nrow(dx), ncol(dx))
  for (j in seq_len(ncol(dx))) for (ids in subsets) {
    boundary <- sort(dx[ids, j])[k]
    closer <- ids[dx[ids, j] < boundary]
    equal <- ids[dx[ids, j] == boundary]
    out[closer, j] <- out[closer, j] + 1 / (length(subsets) * k)
    out[equal, j] <- out[equal, j] + (k - length(closer)) /
      (length(subsets) * k * length(equal))
  }
  out
}

.distance_fixture <- function() list(
  X = matrix(c(-3, -2, -1, 1, 2, 3, 5, 8), ncol = 1L),
  X0 = matrix(c(0, 2, .25), ncol = 1L),
  y = c(3, -2, 8, 1, 7, -5, 4, 9))

.distance_dx <- function(X, X0) vapply(seq_len(nrow(X0)), function(j)
  abs(X[, 1L] - X0[j, 1L]), numeric(nrow(X)))

.distance_capture <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, messages = messages)
}

test_that("prepared group facts are reused without reconstructing geometry", {
  calls <- 0L
  original <- RdistNN:::.dnn_prepare_distance_groups
  local_mocked_bindings(.dnn_prepare_distance_groups = function(dx) {
    calls <<- calls + 1L
    original(dx)
  }, .package = "RdistNN")
  operation <- function() {
    RdistNN:::.dnn_distance_scope()
    prepared <- RdistNN:::.dnn_prepare_distance_groups(matrix(c(1, 1, 2, 3, 4, 5), 3, 2))
    groups <- prepared$groups
    RdistNN:::.dnn_remember_group_facts(groups, prepared$tied_targets)
    local_mocked_bindings(.dnn_compute_group_facts = function(...) stop("group rehash"), .package = "RdistNN")
    expect_identical(RdistNN:::.dnn_group_facts(groups)$tied_targets, c(TRUE, FALSE))
    expect_true(RdistNN:::.dnn_has_ties(groups))
    RdistNN:::.dnn_distance_scope(groups)
    RdistNN:::.dnn_distance_scope(groups)
    RdistNN:::.dnn_average_distance_weights(matrix(1:6, 3, 2), groups)
  }
  for (i in 1:2) expect_warning(operation(), class = "rdistnn_distance_ties_warning")
  expect_identical(calls, 2L)
})

test_that("saved fit validation reuses groups and rejects mutations before notice reuse", {
  d <- .distance_fixture()
  fit <- .with_localization_policy_warnings(np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 3L))
  local_mocked_bindings(.dnn_prepare_distance_groups = function(...) stop("group preparation"), .package = "RdistNN")
  operation <- function() {
    RdistNN:::.dnn_reuse_scope(fit)
    RdistNN:::.dnn_check_report_state(fit)
    expect_true(RdistNN:::.dnn_has_ties(fit$.model$groups))
    broken <- fit
    broken$.model$groups[1, 1] <- broken$.model$groups[3, 1]
    expect_error(RdistNN:::.dnn_check_report_state(broken), class = "rdistnn_state_error")
  }
  expect_warning(operation(), class = "rdistnn_distance_ties_warning")
})

test_that("fresh geometry can become tied within an initially untied operation", {
  operation <- function() {
    RdistNN:::.dnn_distance_scope()
    prepared <- RdistNN:::.dnn_prepare_distance_groups(matrix(1:4))
    first <- prepared$groups
    RdistNN:::.dnn_remember_group_facts(first, prepared$tied_targets)
    expect_false(RdistNN:::.dnn_has_ties(first))
    RdistNN:::.dnn_distance_scope(first)
    prepared <- RdistNN:::.dnn_prepare_distance_groups(matrix(c(1, 1, 3)))
    second <- prepared$groups
    RdistNN:::.dnn_remember_group_facts(second, prepared$tied_targets)
    expect_true(RdistNN:::.dnn_has_ties(second))
    RdistNN:::.dnn_distance_scope(second)
  }
  expect_warning(operation(), class = "rdistnn_distance_ties_warning")
})

test_that("deletion refits do not evict full-sample facts before grid reuse", {
  d <- .distance_fixture()
  scanned <- integer()
  original <- RdistNN:::.dnn_compute_group_facts
  local_mocked_bindings(.dnn_compute_group_facts = function(groups) {
    scanned <<- c(scanned, length(groups))
    original(groups)
  }, .package = "RdistNN")
  .with_localization_policy_warnings(np_dnn(X = d$X, X0 = d$X0, Y = d$y,
    jackknife = "loo", tuning_mode = "grid", return = c("W", "ranks", "dx")))
  expect_false(any(scanned == nrow(d$X) * nrow(d$X0)))
})

test_that("named metric output retains valid group boundaries on serialized reuse", {
  metric <- new_distance("named_distance_output", eval_batch = function(X, X0, ctx = NULL) {
    d <- abs(outer(X[, 1L], X0[, 1L], "-"))
    dimnames(d) <- list(paste0("row", seq_len(nrow(X))), paste0("target", seq_len(nrow(X0))))
    d
  })
  fit <- .with_localization_policy_warnings(np_dnn(X = matrix(1:8), X0 = matrix(c(.5, 2)),
    Y = 1:8, s = 3L, distance = metric, se = FALSE))
  fit <- unserialize(serialize(fit, NULL))
  expect_identical(unname(.with_localization_policy_warnings(coef(fit))), fit$mu)
  expect_identical(.with_localization_policy_warnings(predict(fit,
    X0 = fit$x0[c(2L, 1L, 2L), , drop = FALSE])),
    unname(fit$mu[c(2L, 1L, 2L)]))
})

test_that("candidate weight kernels consume facts without searching for a scope", {
  d <- .distance_fixture()
  original <- RdistNN:::.dnn_grid_compute_key
  local_mocked_bindings(.dnn_grid_compute_key = function(...) {
    local_mocked_bindings(.dnn_reporting_scope = function() stop("scope lookup in candidate kernel"),
      .package = "RdistNN")
    original(...)
  }, .package = "RdistNN")
  fit <- .with_localization_policy_warnings(np_dnn(X = d$X, X0 = d$X0, Y = d$y,
    se = FALSE, tuning_mode = "grid"))
  expect_true(all(fit$tuning_grid$candidates$status == "ok"))
})

test_that("exact distance blocks agree with independent subset enumeration for k1 and general k", {
  d <- .distance_fixture()
  dx <- .distance_dx(d$X, d$X0)
  for (k in c(1L, 2L, 3L)) {
    out <- suppressWarnings(weights_dnn(d$X, d$X0, k = k, s = 5L,
      return = c("ranks", "ord", "dx")))
    expect_equal(out$W, .distance_subset_weights(dx, k, 5L), tolerance = 1e-14)
    expect_identical(as.integer(out$ord[, 1L]), c(3L, 4L, 2L, 5L, 1L, 6L, 7L, 8L))
    expect_identical(dim(out$groups), c(8L, 3L))
    expect_identical(out$groups[3L, 1L], out$groups[4L, 1L])
    expect_false(out$groups[3L, 3L] == out$groups[4L, 3L])
  }
})

test_that("distance metadata is preserved with authoritative supplied ranks and bare ranks stay strict", {
  X <- matrix(c(-1, 1, -1, 2, 4), ncol = 1L)
  ranks <- matrix(c(5L, 1L, 3L, 2L, 4L), ncol = 1L)
  dx <- matrix(c(1, 1, 1, 2, 4), ncol = 1L)
  rw <- as_rank_weight(function(r) c(2, -3, 5, 7, -11)[r])
  raw <- matrix(c(-11, 2, 5, -3, 7), ncol = 1L)
  strict <- weights_matrix(X, x0 = 0, ranking = ranks, weight = rw)
  expect_identical(strict$W, raw)
  for (ranking in list(list(ranks = ranks, dx = dx),
      list(ord = matrix(order(ranks[, 1L]), ncol = 1L), dx = dx))) {
    out <- suppressWarnings(weights_matrix(X, x0 = 0, ranking = ranking,
      weight = rw, return = c("ranks", "dx")))
    expected <- raw
    expected[1:3, 1L] <- sum(raw[1:3, 1L]) / 3
    expect_identical(out$ranks, ranks)
    expect_identical(out$dx, dx)
    expect_equal(out$W, expected, tolerance = 0)
    expect_equal(sum(out$W), sum(raw), tolerance = 1e-15)
  }
})

test_that("exact equality uses represented distances without a near-tie tolerance", {
  X <- matrix(c(-1, 1, 1 + 2^-51, 2), ncol = 1L)
  out <- suppressWarnings(weights_dnn(X, x0 = 0, k = 1L, s = 3L, return = "dx"))
  expected <- .distance_subset_weights(matrix(abs(X[, 1L]), ncol = 1L), 1L, 3L)
  expect_equal(out$W, expected, tolerance = 1e-15)
  expect_equal(out$W[1L], out$W[2L])
  expect_gt(out$W[2L], out$W[3L])
})

test_that("distance block fast and direct deletion engines retain actual survivor rows", {
  d <- .distance_fixture()
  dx <- .distance_dx(d$X, d$X0)
  fit <- suppressWarnings(weights_dnn(d$X, d$X0, k = 2L, s = 5L, return = "ranks"))
  signal <- cbind(d$y, rev(d$y))
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  expect_true(is.function(RdistNN:::.rdistnn_compiled_fun("rustat_block_rank_loo")))
  reference <- array(0, c(8L, 2L, 3L))
  for (i in 1:8) reference[i, , ] <- crossprod(signal[-i, , drop = FALSE],
    .distance_subset_weights(dx[-i, , drop = FALSE], 2L, 5L))
  for (compiled in c(TRUE, FALSE)) for (direct in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    rw <- rank_weight_dnn_exact(2L, 5L)
    jk <- RdistNN:::.dnn_loo_centered(signal, fit$ranks,
      .block_exact_table(7L, 2L, 5L), .block_exact_table(8L, 2L, 5L),
      rw = if (direct) rw else NULL, groups = fit$groups)
    expect_equal(RdistNN:::.dnn_restore_loo(jk$loo, jk$offset), reference, tolerance = 1e-12)
  }
})

test_that("target distance groups survive serialization, reflection and candidate reconstruction", {
  d <- .distance_fixture()
  tuning <- new_tuning("distance013", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 5L)))
  f <- suppressWarnings(np_dnn(X = d$X, X0 = d$X0, Y = d$y,
    tuning = tuning, tuning_mode = "grid", return = c("W", "ranks", "ess")))
  expect_identical(dim(f$.model$groups), c(8L, 3L))
  reflected <- suppressWarnings(np_dnn(X = -d$X, X0 = -d$X0, Y = d$y,
    tuning = tuning, tuning_mode = "grid", return = c("W", "ranks", "ess")))
  expect_equal(f$W, reflected$W, tolerance = 1e-14)
  expect_equal(f$mu, reflected$mu, tolerance = 1e-13)
  expect_equal(f$se, reflected$se, tolerance = 1e-13)
  f <- unserialize(serialize(f, NULL))
  for (i in 1:3) {
    view <- suppressWarnings(candidate_fit(f, i))
    W <- .distance_subset_weights(.distance_dx(d$X, d$X0), 1L, view$meta$s)
    expect_equal(view$W, W, tolerance = 1e-14)
    expect_equal(view$mu, as.numeric(crossprod(d$y, W)), tolerance = 1e-13)
    expect_equal(suppressWarnings(diagnostics(view))$n_eff, 1 / colSums(W^2), tolerance = 1e-13)
  }
  expect_equal(suppressWarnings(predict(f, X0 = d$X0, se.fit = TRUE))$fit,
    f$mu, tolerance = 1e-13)
  invalid_policy <- f
  invalid_policy$.model$tie_policy <- "distance_then_covariates_with_duplicate_block_averaging"
  expect_error(predict(invalid_policy, X0 = d$X0), "Refit|refit")
})

test_that("public warning reports actual distance ties once even discovered after an untied batch", {
  d <- .distance_fixture()
  quiet <- .distance_capture(weights_dnn(d$X, x0 = .25, k = 1L, s = 5L))
  expect_length(quiet$messages, 0L)
  warned <- .distance_capture(np_dnn(X = d$X, X0 = d$X0, Y = d$y,
    s = 5L, jackknife = "loo"))
  tie <- warned$messages[grepl("continuous distance distribution", warned$messages, fixed = TRUE)]
  expect_length(tie, 1L)
  expect_match(tie, "persistent ties", fixed = TRUE)
  outer <- function() {
    RdistNN:::.dnn_distance_scope()
    weights_dnn(d$X, x0 = .25, k = 1L, s = 5L)
    weights_dnn(d$X, x0 = 0, k = 1L, s = 5L)
  }
  expect_length(.distance_capture(outer())$messages, 1L)
})

test_that("saved target groups are validated against retained represented distances", {
  d <- .distance_fixture()
  f <- suppressWarnings(np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 5L,
    return = character()))
  expect_identical(dim(f$.model$ranks), c(8L, 3L))
  expect_identical(dim(f$.model$dx), c(8L, 3L))
  corrupt <- f
  corrupt$.model$groups[, 1L] <- seq_len(8L)
  expect_error(coef(corrupt), "distance groups|Refit|refit")
})

test_that("malformed saved distance group state fails cleanly before coercion or member access", {
  d <- .distance_fixture()
  f <- suppressWarnings(np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 5L,
    return = character()))
  for (bad in c(2^32, Inf, NA_real_, -1, .5)) {
    corrupt <- f
    corrupt$.model$groups[1L, 1L] <- bad
    expect_warning(expect_error(coef(corrupt), "distance groups|Refit|refit"), NA)
  }
  corrupt <- f
  corrupt$.model$ranks <- corrupt$.model$dx <- NULL
  expect_error(coef(corrupt), "distance groups|Refit|refit")
})

test_that("all scalar adapters retain one computed group cache without rerunning distance callbacks", {
  d <- .distance_fixture()
  calls <- 0L
  locked <- FALSE
  metric <- new_distance("distance013_stateful", eval_batch = function(X, X0, ctx) {
    if (locked) stop("distance callback must not rerun for retained targets")
    calls <<- calls + 1L
    .distance_dx(X, X0)
  })
  tuning <- new_tuning("distance013_cache", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 5L)))
  make <- function(kind) {
    common <- list(X = d$X, X0 = d$X0, distance = metric,
      tuning = tuning, tuning_mode = "grid", return = c("W", "ess"))
    switch(kind,
      np = do.call(np_dnn, c(common, list(Y = d$y))),
      cate = do.call(cate_dnn, c(common, list(pseudo_outcome = d$y))),
      ratio = do.call(affine_ratio_dnn, c(common, list(numerator = d$y,
        denominator = seq(.5, 1.2, length.out = 8L)))),
      plqr = do.call(plqr_dnn, c(common, list(Y = d$y,
        D = matrix(seq_len(8L)), mD = matrix(0, 8L), g = rep(0, 8L),
        theta_eval = .2, tau = .5))))
  }
  for (kind in c("np", "cate", "ratio", "plqr")) {
    locked <- FALSE
    before <- calls
    f <- suppressWarnings(make(kind))
    expect_identical(calls - before, 1L, info = kind)
    expect_identical(dim(f$.model$ranks), c(8L, 3L))
    expect_identical(dim(f$.model$groups), c(8L, 3L))
    locked <- TRUE
    restored <- unserialize(serialize(f, NULL))
    views <- suppressWarnings(candidate_fits(restored))
    for (view in views) {
      expect_equal(suppressWarnings(diagnostics(view))$n_eff, 1 / colSums(view$W^2))
      expect_s3_class(suppressWarnings(plot_data_weights(view)), "data.frame")
      if (kind != "plqr") expect_true(all(is.finite(unlist(suppressWarnings(vcov(view))))))
    }
    expect_identical(calls - before, 1L, info = kind)
    locked <- FALSE
    suppressWarnings(predict(f, X0 = d$X0[3:1, , drop = FALSE], se.fit = kind != "plqr"))
    expect_identical(calls - before, 2L, info = kind)
  }
})

test_that("deletion-refitted geometry rebuilds distance groups while fast and ratio paths stay frozen", {
  d <- .distance_fixture()
  metric <- new_distance("distance013_refit_groups",
    setup = function(X) list(reduced = nrow(X) < 8L),
    eval_batch = function(X, X0, ctx) {
      dx <- .distance_dx(X, X0)
      if (ctx$reduced) floor(dx / 2) else dx
    })
  dx <- .distance_dx(d$X, d$X0)
  D <- matrix(seq_len(8L))
  score <- D * (.5 - as.numeric(d$y <= .2 * D[, 1L]))
  deletion <- function(signal, refit, tdnn = FALSE) {
    t(vapply(seq_len(8L), function(i) {
      reduced_dx <- dx[-i, , drop = FALSE]
      if (refit) reduced_dx <- floor(reduced_dx / 2)
      W <- .distance_subset_weights(reduced_dx, 1L, 3L)
      if (tdnn) {
        b <- (5 / 2)^2
        W <- (b * .distance_subset_weights(reduced_dx, 1L, 5L) -
          .distance_subset_weights(reduced_dx, 1L, 2L)) / (b - 1)
      }
      as.numeric(crossprod(signal[-i], W))
    }, numeric(3L)))
  }
  for (engine in c("fast", "loo")) for (kind in c("np", "cate", "tdnn", "plqr")) {
    common <- list(X = d$X, X0 = d$X0, distance = metric, jackknife = engine,
      return = character())
    f <- suppressWarnings(switch(kind,
      np = do.call(np_dnn, c(common, list(Y = d$y, s = 3L))),
      cate = do.call(cate_dnn, c(common, list(pseudo_outcome = d$y, s = 3L))),
      tdnn = do.call(np_tdnn, c(common, list(Y = d$y, s1 = 2L, s2 = 5L, dim = 1L))),
      plqr = do.call(plqr_dnn, c(common, list(Y = d$y, D = D,
        mD = matrix(0, 8L), g = rep(0, 8L), theta_eval = .2, tau = .5, s = 3L)))))
    L <- deletion(if (kind == "plqr") as.numeric(score) else d$y,
      refit = engine == "loo", tdnn = kind == "tdnn")
    expected <- 7 / 8 * colSums(sweep(L, 2L, colMeans(L), "-")^2)
    actual <- if (kind == "plqr") as.numeric(f$moment_se)^2 else f$se^2
    expect_equal(actual, expected, tolerance = 2e-12, info = paste(kind, engine))
  }
  # The paired ratio explicit route freezes the same original distance groups.
  for (engine in c("fast", "loo")) {
    f <- suppressWarnings(affine_ratio_dnn(X = d$X, X0 = d$X0,
      numerator = d$y, denominator = rev(d$y) + 10, distance = metric,
      s = 3L, jackknife = engine, return = "deletions"))
    expect_equal(f$deletions[, 1L, ], deletion(d$y, FALSE), tolerance = 1e-12)
    expect_equal(f$deletions[, 2L, ], deletion(rev(d$y) + 10, FALSE), tolerance = 1e-12)
  }
  # The full sample is untied here; ties first arise inside deletion geometry.
  capture <- .distance_capture(np_dnn(X = d$X, x0 = .25, Y = d$y,
    distance = metric, s = 3L, jackknife = "loo"))
  expect_length(capture$messages[grepl("continuous distance distribution", capture$messages, fixed = TRUE)], 1L)
})

test_that("GMM averages target distance weights before solving and retains the same geometry", {
  d <- .distance_fixture()
  signal <- cbind(Y1 = d$y, Y2 = rev(d$y)^2)
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  expect_true(is.function(RdistNN:::.rdistnn_compiled_fun("rustat_block_rank_loo")))
  for (compiled in c(TRUE, FALSE)) for (k in c(1L, 2L)) {
    options(RdistNN.use_compiled = compiled)
    W <- .distance_subset_weights(.distance_dx(d$X, d$X0), k, 5L)
    fit <- suppressWarnings(gmm_dnn(X = d$X, X0 = d$X0,
      model_data = as.data.frame(signal), model = .test_mean2_gmm_model(),
      k = k, s = 5L, return = c("W", "ess")))
    expect_equal(fit$W, W, tolerance = 1e-14)
    expect_equal(unname(fit$theta), unname(crossprod(signal, W)), tolerance = 1e-10)
    expect_identical(dim(fit$.model$groups), c(8L, 3L))
    expect_identical(dim(fit$.model$ranks), c(8L, 3L))
    saved <- unserialize(serialize(fit, NULL))
    expect_equal(unname(saved$neff), 1 / colSums(W^2))
    expect_equal(suppressWarnings(plot_data_weights(saved))$weight,
      as.numeric(vapply(seq_len(3L), function(j) W[order(fit$.model$ranks[, j]), j], numeric(8L))))
  }
})

test_that("nested custom ranking callbacks share the public weight warning scope", {
  d <- .distance_fixture()
  ranker <- function(X, X0, ctx) {
    weights_dnn(X, X0, k = 1L, s = 5L)
    list(dx = .distance_dx(X, X0))
  }
  captured <- .distance_capture(weights_matrix(d$X, d$X0, ranking = ranker,
    weight = rank_weight_dnn_exact(1L, 5L)))
  expect_length(captured$messages, 1L)
})

test_that("the TDNN fast adapter calls the linear block kernel", {
  d <- .distance_fixture()
  wt <- suppressWarnings(weights_tdnn(d$X, d$X0, s1 = 2L, s2 = 5L,
    dim = 1L, return = "ranks"))
  paths <- logical()
  real_block_loo <- RdistNN:::.dnn_block_loo
  testthat::local_mocked_bindings(.dnn_block_loo = function(signal, ranks, table,
      groups, direct = FALSE) {
    paths <<- c(paths, direct)
    real_block_loo(signal, ranks, table, groups, direct)
  }, .package = "RdistNN")
  RdistNN:::.tdnn_jackknife_se_loo_matrix(d$X, d$X0, d$y, "l2", "exact",
    2L, 5L, 1L, TRUE, NULL, jackknife = "fast",
    precomputed_ranks = wt$ranks, precomputed_groups = wt$groups)
  expect_identical(paths, c(FALSE, FALSE))
})

test_that("current saved fits require canonical finite rank permutations and complete caches", {
  d <- .distance_fixture()
  fit <- suppressWarnings(np_dnn(X = d$X, X0 = d$X0, Y = d$y, s = 5L,
    return = character()))
  corruptions <- list(
    reversed = function(model) { model$ranks <- model$ranks[8:1, , drop = FALSE]; model },
    tied_swap = function(model) { model$ranks[c(3L, 4L), 1L] <- model$ranks[c(4L, 3L), 1L]; model },
    missing = function(model) { model$ranks <- NULL; model },
    duplicated = function(model) { model$ranks[1L, 1L] <- model$ranks[2L, 1L]; model },
    nonfinite = function(model) { model$ranks[1L, 1L] <- NA_integer_; model },
    overflow = function(model) { model$ranks[1L, 1L] <- 2^32; model },
    missing_x = function(model) { model$x_mat <- NULL; model },
    missing_dx = function(model) { model$dx <- NULL; model },
    missing_model = function(model) NULL)
  for (name in names(corruptions)) {
    bad <- fit
    bad$.model <- corruptions[[name]](bad$.model)
    for (consumer in list(coef, vcov, summary, plot_data_weights,
                         function(x) predict(x, x0 = .25))) {
      expect_warning(expect_error(consumer(bad), "Refit|refit", info = name), NA)
    }
  }
})
