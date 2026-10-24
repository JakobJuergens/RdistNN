test_that("CLATE candidate views preserve paired inference and frozen state", {
  for (method in c("exact", "asymp")) for (normalize in c(TRUE, FALSE)) {
    a <- .clate_fit_arguments()
    a$s <- NULL
    a$tuning <- .clate_alignment_tuning()
    a$return <- c("W", "ess", "neighbors", "deletions", "pseudo_values")
    expect_warning(grid <- do.call(clate_dnn, c(a, list(method = method,
      k = 2L, normalize = normalize, tuning_mode = "grid"))), "infeasible")
    candidates <- grid$tuning_grid$candidates
    ids <- candidates$candidate_id[candidates$feasible]
    views <- candidate_fits(grid, ids)
    expect_named(views, ids)
    for (i in which(candidates$feasible)) {
      args <- a
      args$tuning <- NULL
      args$s <- candidates$s[i]
      scalar <- suppressWarnings(do.call(clate_dnn,
        c(args, list(method = method, k = 2L, normalize = normalize))))
      view <- candidate_fit(grid, i)
      expect_s3_class(view, "rdistnn_clate_result")
      expect_equal(view, views[[candidates$candidate_id[i]]])
      for (field in c("theta", "se", "components", "component_vcov", "theta_vcov",
                      "theta_omega", "status", "neff", "W", "ess", "neighbors",
                      "deletions", "pseudo_values")) {
        expect_equal(view[[field]], scalar[[field]], tolerance = 1e-12, info = field)
      }
      expect_identical(view$.model$iv, grid$.model$iv)
      expect_equal(summary(view)$coefficients, summary(scalar)$coefficients)
      expect_equal(suppressWarnings(predict(view, X0 = a$X0, se.fit = TRUE)),
                   suppressWarnings(predict(scalar, X0 = a$X0, se.fit = TRUE)))
      expect_null(view$tuning_grid)
    }
    expect_error(candidate_fit(grid, which(!candidates$feasible)), "infeasible")
  }
})

test_that("point-only ratio candidate views work without public extras or ranker calls", {
  calls <- 0L
  distance <- new_distance("candidate_frozen", eval = function(X, x0, ctx = NULL) {
    calls <<- calls + 1L
    abs(X[, 1L] - x0[1L])
  })
  fit <- affine_ratio_dnn(X = matrix(1:8), x0 = 0, numerator = 1:8,
    denominator = rep(1, 8), tuning = .clate_alignment_tuning(),
    tuning_mode = "grid", se = FALSE, return = character(), distance = distance)
  view <- candidate_fit(fit, 1L)
  expect_equal(calls, 1L)
  expect_null(view$se)
  expect_null(view$W)
  # Uniformly average the minima of all 28 size-two subsets of 1:8.
  expect_equal(view$theta, 3)
  expect_false(any(view$status$wald_available))
})
