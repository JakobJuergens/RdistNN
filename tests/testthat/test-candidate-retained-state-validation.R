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
    suppressWarnings(do.call(fn,
      c(common, list(tuning = tuning, tuning_mode = "grid"))))
  })
  list(fits = fits, constructors = constructors, common = common)
}

test_that("candidate payloads require nonempty unique candidate IDs", {
  fits <- .candidate_retained_state_fixtures()$fits
  for (fit in fits) {
    for (invalid in list(NULL, rep(NA_character_, 4L), rep("", 4L), rep("same", 4L))) {
      bad <- fit
      bad$tuning_grid$candidates$candidate_id <- invalid
      expect_error(suppressWarnings(candidate_fit(bad, 1L)), "candidate_id.*(non-empty|unique)")
    }
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
      expect_error(suppressWarnings(candidate_fit(bad, 1L)), "integer.*target_id")
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
    unavailable$tuning_grid$candidates$reason[4L] <- "unavailable order"
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
      if (!is.null(bad$.model$ranks)) {
        bad$.model$ranks <- invalid
      } else {
        bad$.model$geometry$ranks <- invalid
      }
      expect_error(suppressWarnings(candidate_fit(bad, 1L)), "[Rr]anks|permutation")
    }
  }
})

test_that("serialized thin candidate views preserve valid duplicate-order rows", {
  fixture <- .candidate_retained_state_fixtures()
  value <- function(fit) {
    if (inherits(fit, "plqr_result")) fit$moment else stats::coef(fit)
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
