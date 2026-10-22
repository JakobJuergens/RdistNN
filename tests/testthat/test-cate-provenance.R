.cate_provenance_case <- function() {
  x <- seq(-1, 1, length.out = 12L)
  data.frame(x = x, y = x^2 + rep(c(-.2, .4), 6L),
    d = rep(0:1, 6L), e = seq(.2, .8, length.out = 12L),
    m0 = x^2, m1 = x^2 + .6)
}

.cate_provenance_fit <- function(d, ...) {
  cate_dnn(data = d, x_cols = "x", y_col = "y", t_col = "d",
    e_col = "e", m0_col = "m0", m1_col = "m1", ...)
}

test_that("CATE source records the input branch rather than signal values or names", {
  d <- .cate_provenance_case()
  d$aipw_oracle <- with(d, m1 - m0 + d * (y - m1) / e -
    (1 - d) * (y - m0) / (1 - e))
  version <- as.character(getNamespaceVersion("RdistNN"))
  for (targets in list(matrix(0), matrix(c(-.3, .4), ncol = 1L))) {
    built <- .cate_provenance_fit(d, X0 = targets, s = 3L, return = "psi")
    matrix_fit <- cate_dnn(X = as.matrix(d["x"]), X0 = targets,
      pseudo_outcome = d$aipw_oracle, s = 3L)
    supplied <- cate_dnn(data = d, x_cols = "x", X0 = targets,
      pseudo_outcome = d$aipw_oracle, s = 3L)
    precedence <- .cate_provenance_fit(d, X0 = targets, s = 3L,
      pseudo_outcome = d$aipw_oracle, return = "psi")
    expect_identical(built$meta$signal_source, "aipw")
    expect_equal(built$psi, d$aipw_oracle, tolerance = 1e-14)
    for (fit in list(matrix_fit, supplied, precedence)) {
      expect_identical(fit$meta$signal_source, "supplied")
      expect_identical(fit$meta$created_with_version, version)
      expect_identical(fit$tau, built$tau)
      expect_identical(fit$se, built$se)
      expect_identical(vcov(fit), vcov(built))
    }
    expect_identical(built$meta$created_with_version, version)
    expect_identical(precedence$psi, d$aipw_oracle)
  }
  ignored <- cate_dnn(data = d, x_cols = "x", x0 = 0, s = 3L,
    pseudo_outcome = d$aipw_oracle, m0_col = "absent")
  expect_identical(ignored$meta$signal_source, "supplied")
})

test_that("CATE summaries and restored fits preserve descriptive metadata", {
  d <- .cate_provenance_case()
  fit <- .cate_provenance_fit(d, x0 = 0, s = 3L)
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(restored$meta, fit$meta)
  expect_identical(summary(restored)$meta, fit$meta)
  expect_output(print(restored), "signal.*aipw")
  expect_output(print(summary(restored)), "signal.*aipw")
  supplied <- cate_dnn(X = as.matrix(d["x"]), x0 = 0,
    pseudo_outcome = d$y, s = 3L)
  expect_output(print(supplied), "signal.*supplied")
  expect_output(print(summary(supplied)), "signal.*supplied")
  expect_identical(coef(restored), coef(fit))
  expect_identical(vcov(restored), vcov(fit))
  expect_identical(predict(restored, x0 = .2, se.fit = TRUE),
    predict(fit, x0 = .2, se.fit = TRUE))
})

test_that("CATE candidates retain original signal source and creation version", {
  d <- .cate_provenance_case()
  tuning <- new_tuning("provenance", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 4L)))
  for (source in c("aipw", "supplied")) {
    fit <- .cate_provenance_fit(d, X0 = matrix(c(-.3, .4), ncol = 1L),
      tuning = tuning, tuning_mode = "grid",
      pseudo_outcome = if (source == "supplied") d$y else NULL)
    expect_identical(fit$meta$signal_source, source)
    fit$meta$created_with_version <- "0.1.0.1"
    for (i in seq_along(candidate_fits(fit))) {
      candidate <- candidate_fit(fit, i)
      expect_identical(candidate$meta$signal_source, source)
      expect_identical(candidate$meta$created_with_version, "0.1.0.1")
      expect_identical(summary(candidate)$meta$created_with_version, "0.1.0.1")
    }
  }
})

test_that("current CATE fits preserve prediction shapes and serialized state", {
  d <- .cate_provenance_case()
  for (jackknife in c("fast", "loo")) {
    fit <- .cate_provenance_fit(d, X0 = matrix(c(-.3, .4), ncol = 1L),
      s = 3L, jackknife = jackknife, distance = distance_l2_scaled("sd"),
      return = c("psi", "nuisance", "W", "ranks"))
    restored <- unserialize(serialize(fit, NULL))
    before <- restored
    expect_identical(coef(restored), coef(fit))
    expect_identical(restored$se, fit$se)
    expect_identical(vcov(restored), vcov(fit))
    expect_identical(confint(restored), confint(fit))
    expect_identical(summary(restored)$coefficients, summary(fit)$coefficients)
    for (m in c(1L, 2L)) {
      targets <- matrix(c(-.1, .2)[seq_len(m)], ncol = 1L)
      for (se.fit in c(FALSE, TRUE)) for (interval in c("none", "confidence")) {
        prediction <- predict(restored, X0 = targets, se.fit = se.fit, interval = interval)
        expect_identical(prediction,
          predict(fit, X0 = targets, se.fit = se.fit, interval = interval))
        if (se.fit) {
          expect_identical(names(prediction), c("fit", "se.fit"))
          expect_identical(attributes(prediction), list(names = c("fit", "se.fit")))
          expect_null(attributes(prediction$se.fit))
          prediction <- prediction$fit
        }
        if (interval == "none") {
          expect_type(prediction, "double")
          expect_length(prediction, m)
          expect_null(attributes(prediction))
        } else {
          expect_identical(dim(prediction), c(m, 3L))
          expect_identical(names(attributes(prediction)), c("dim", "dimnames"))
        }
      }
    }
    .expect_identical_fast(restored, before)
    restored$.model <- NULL
    expect_error(predict(restored, x0 = 0), "fitted model state")
    expect_error(vcov(restored), "fitted model state")
  }
})

test_that("retired pooled-mean CATE arguments remain rejected", {
  d <- .cate_provenance_case()
  expect_error(cate_dnn(data = d, x_cols = "x", y_col = "y", t_col = "d",
    e_col = "e", m_col = "m0", x0 = 0, s = 3L), "unused argument.*m_col")
})

test_that("CATE reporting requires current provenance metadata", {
  fit <- .cate_provenance_fit(.cate_provenance_case(), x0 = 0, s = 3L)
  for (field in c("signal_source", "created_with_version")) {
    bad <- fit
    bad$meta[[field]] <- NULL
    expect_error(summary(bad), "metadata.*[Rr]efit")
    expect_error(print(bad), "metadata.*[Rr]efit")
  }
})
