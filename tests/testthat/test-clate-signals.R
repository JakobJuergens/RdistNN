test_that("augmented responses match explicit formulas and share propensity", {
  a <- .clate_test_predictions()
  z <- do.call(iv_augmented_signals, a)
  want_y <- a$m_y1 - a$m_y0 + a$A / a$e * (a$Y - a$m_y1) -
    (1 - a$A) / (1 - a$e) * (a$Y - a$m_y0)
  want_d <- a$m_d1 - a$m_d0 + a$A / a$e * (a$D - a$m_d1) -
    (1 - a$A) / (1 - a$e) * (a$D - a$m_d0)
  expect_s3_class(z, "iv_signals")
  expect_equal(z$phi_y, want_y)
  expect_equal(z$phi_d, want_d)
  expect_identical(z$row_ids, seq_len(8L))
  expect_identical(z$propensity$raw, a$e)
  expect_identical(z$propensity$used, a$e)
  expect_false(any(z$propensity$clipped))
})

test_that("one-sided and unrestricted signals agree without censoring negatives", {
  a <- .clate_test_predictions()
  a$D <- a$D * a$A
  a$m_d0 <- rep(0, 8L)
  general <- do.call(iv_augmented_signals, a)
  a$m_d0 <- NULL
  restricted <- do.call(iv_augmented_signals, c(a, list(one_sided = TRUE)))
  expect_equal(restricted$phi_d, general$phi_d)
  expect_equal(restricted$phi_y, general$phi_y)
  expect_equal(restricted$nuisance$m_d0, rep(0, 8L))
  z <- iv_augmented_signals(
    Y = 0, D = 0, A = 1, e = .5,
    m_y0 = 0, m_y1 = 0, m_d1 = .6, one_sided = TRUE
  )
  expect_equal(z$phi_d, -.6)
})

test_that("perfect compliance and full instrument relabeling preserve identities", {
  a <- .clate_test_predictions()
  a$D <- a$A
  a$m_d0 <- rep(0, 8L)
  a$m_d1 <- rep(1, 8L)
  z <- do.call(iv_augmented_signals, a)
  expect_equal(z$phi_d, rep(1, 8L))
  flipped <- a
  flipped$A <- 1 - a$A
  flipped$e <- 1 - a$e
  flipped$m_y0 <- a$m_y1
  flipped$m_y1 <- a$m_y0
  flipped$m_d0 <- a$m_d1
  flipped$m_d1 <- a$m_d0
  reverse <- do.call(iv_augmented_signals, flipped)
  expect_equal(reverse$phi_y, -z$phi_y)
  expect_equal(reverse$phi_d, -z$phi_d)
})

test_that("outcome location and scale transformations affect only the contrast", {
  a <- .clate_test_predictions()
  original <- do.call(iv_augmented_signals, a)
  for (scale in c(-3, 0, 2)) {
    b <- a
    for (nm in c("Y", "m_y0", "m_y1")) b[[nm]] <- 7 + scale * a[[nm]]
    transformed <- do.call(iv_augmented_signals, b)
    expect_equal(transformed$phi_y, scale * original$phi_y, tolerance = 1e-12)
    expect_equal(transformed$phi_d, original$phi_d)
  }
})

test_that("propensity endpoints are clipped once and provenance is retained", {
  a <- .clate_test_predictions()
  a$e[c(1L, 2L)] <- c(0, 1)
  z <- do.call(iv_augmented_signals, c(a, list(clip = .1)))
  expect_equal(z$propensity$used, c(.1, .9, a$e[-c(1L, 2L)]))
  expect_identical(z$propensity$raw, a$e)
  expect_equal(which(z$propensity$clipped), c(1L, 2L))
  expect_true(all(is.finite(c(z$phi_y, z$phi_d))))
})

test_that("product remainder has both positive regression-error terms", {
  for (e0 in c(.2, .5, .8)) for (e in c(.3, .6)) {
    truth <- c(-.7, 1.3)
    for (error in list(c(.4, -.2), c(0, 0), c(-.3, .5))) {
      fit <- truth + error
      # Enumerate assignment at its true probability; deterministic arm responses.
      s <- iv_augmented_signals(
        Y = truth, D = c(0, 1), A = c(0, 1), e = rep(e, 2L),
        m_y0 = rep(fit[1L], 2L), m_y1 = rep(fit[2L], 2L),
        m_d0 = rep(0, 2L), m_d1 = rep(1, 2L)
      )
      drift <- sum(c(1 - e0, e0) * s$phi_y) - diff(truth)
      want <- (e - e0) * (error[2L] / e + error[1L] / (1 - e))
      expect_equal(as.numeric(drift), as.numeric(want), tolerance = 1e-12)
    }
  }
})

test_that("signal construction rejects malformed vectors and invalid semantics", {
  a <- .clate_test_predictions()
  for (nm in names(a)) {
    for (bad in list(NA_real_, NaN, Inf, -Inf)) {
      b <- a
      b[[nm]][2L] <- bad
      expect_error(do.call(iv_augmented_signals, b), nm)
    }
    b <- a
    b[[nm]] <- b[[nm]][-1L]
    expect_error(do.call(iv_augmented_signals, b), if (nm == "Y") "D" else nm)
    b[[nm]] <- matrix(a[[nm]], ncol = 1L)
    expect_error(do.call(iv_augmented_signals, b), nm)
  }
  for (nm in c("D", "A")) for (bad in list(rep(TRUE, 8L), factor(a[[nm]]), rep(2, 8L))) {
    b <- a
    b[[nm]] <- bad
    expect_error(do.call(iv_augmented_signals, b), nm)
  }
  for (nm in c("e", "m_d0", "m_d1")) {
    b <- a
    b[[nm]][1L] <- -0.01
    expect_error(do.call(iv_augmented_signals, b), nm)
  }
  for (clip in list(0, .5, NA_real_, c(.1, .2), "x")) {
    expect_error(do.call(iv_augmented_signals, c(a, list(clip = clip))), "clip")
  }
  expect_error(do.call(iv_augmented_signals, c(a, list(one_sided = TRUE))), "one_sided")
  a$m_d0 <- NULL
  expect_error(do.call(iv_augmented_signals, a), "m_d0")
})

test_that("row identities and training declarations are explicit", {
  a <- .clate_test_predictions()
  ids <- paste0("household_", seq_len(8L))
  z <- do.call(iv_augmented_signals, c(a, list(row_ids = ids, training = "full_sample")))
  expect_identical(z$row_ids, ids)
  expect_identical(z$meta$training, "full_sample")
  for (bad in list(rep("a", 8L), c(ids[-1L], NA_character_),
                   c(ids[-1L], ""), 1:7, c(1:7, Inf), factor(ids))) {
    expect_error(do.call(iv_augmented_signals, c(a, list(row_ids = bad))), "row_ids")
  }
  expect_error(do.call(iv_augmented_signals, c(a, list(training = "auto"))), "training")
})
