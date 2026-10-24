test_that("ratio diagnostics reconstruct paired uncertainty and score concentration", {
  f <- .ratio_fixture()
  z <- diagnostics(f)
  expect_s3_class(z, "rdistnn_diagnostics")
  expect_equal(z$numerator, 2)
  expect_equal(z$denominator, .5)
  expect_equal(z$numerator_se, sqrt(.5))
  expect_equal(z$denominator_se, sqrt(.75))
  expect_equal(z$component_correlation, .5 / sqrt(.5 * .75))
  expect_equal(z$component_covariance, .5)
  expect_identical(z$geometry_rule, "frozen_full_sample")
  expect_identical(z$tie_rule, "exact_distance_block_averaging_v1")
  expect_equal(z$neff, 18 / 7)
  expect_equal(z$max_weight, .5)
  expect_equal(z$radius_90, 3)
  expect_equal(z$jk_top1_sq_share, 64 / 102)
  expect_equal(z$jk_top5_sq_share, 1)
  expect_false(z$instrument_available)
  expect_true(is.na(z$instrument_1_mass))
  expect_true(is.na(z$propensity_clipped_mass))
  expect_error(diagnostics(f, candidates = "all"), "grid|tuning_mode")
})

test_that("CLATE diagnostics preserve IV roles clipping and instrument support", {
  a <- .clate_fit_arguments()
  a$e[1L] <- 0
  a$return <- c("W", "pseudo_values")
  f <- do.call(clate_dnn, c(a, list(clip = .1, training = "oob")))
  z <- diagnostics(f)
  expect_true(all(z$instrument_available))
  expect_equal(z$instrument_1_mass, as.numeric(crossprod(a$A, f$W)))
  expect_equal(z$instrument_0_mass, 1 - z$instrument_1_mass)
  expect_equal(z$propensity_clipped_mass, as.numeric(f$W[1L, ]))
  expect_identical(z$training, rep("oob", 2L))
  expect_identical(z$numerator_role, rep("reduced_form", 2L))
  expect_identical(z$denominator_role, rep("first_stage", 2L))
  pv <- f$pseudo_values
  for (j in 1:2) {
    q <- f$theta[j] * pv[, 2L, j] - pv[, 1L, j]
    h <- (q - mean(q))^2
    expect_equal(z$jk_top1_sq_share[j], max(h) / sum(h), tolerance = 1e-12)
  }
})

test_that("diagnostics handle deterministic denominators and disabled covariance", {
  f <- .ratio_fixture(denominator = rep(1, 4L))
  z <- diagnostics(f)
  expect_equal(z$denominator_se, 0)
  expect_true(is.na(z$denominator_stat))
  expect_identical(z$denominator_variance_status, "zero")
  expect_true(is.na(z$component_correlation))
  point <- diagnostics(.ratio_fixture(se = FALSE))
  expect_true(is.na(point$se))
  expect_true(is.na(point$jk_top1_sq_share))
  expect_false(point$score_concentration_available)
  expect_true(is.na(point$component_covariance))
})
test_that("paired diagnostics work without utils attached to the search path", {
  fit <- .ratio_fixture()
  without_utils <- function(code) {
    position <- match("package:utils", search())
    if (!is.na(position)) {
      detach("package:utils", character.only = TRUE)
      on.exit(suppressPackageStartupMessages(library("utils", character.only = TRUE, pos = position)), add = TRUE)
    }
    force(code)
  }
  expect_s3_class(without_utils(diagnostics(fit)), "rdistnn_diagnostics")
})
