.retained_field_fixture <- function(family, frame = FALSE, supplied = FALSE,
                                    x_cols = c("a", "b")) {
  a <- seq_len(12L)
  X <- cbind(a = a, b = a^2 / 10)
  arm <- rep(c(0, 1), 6L)
  d <- data.frame(a = a, b = X[, 2L], y = a^2, observed = arm * a^2,
    arm = arm, e = seq(.15, .85, length.out = 12L), m0 = 0, m1 = a^2,
    d0 = 0, d1 = 1, den = 1)
  d$pair <- I(X)
  common <- list(X0 = rbind(c(.2, 5), c(10, .5)), method = "asymp", se = TRUE,
    return = c("W", "ranks", "dx"))
  common <- c(common, if (family == "tdnn") list(s1 = 2L, s2 = 4L, dim = 2L) else list(k = 2L, s = 4L))
  args <- if (frame) c(list(data = d, x_cols = x_cols), switch(family,
    np = list(y_col = "y"), tdnn = list(y_col = "y"),
    cate = if (supplied) list(pseudo_outcome = d$y) else list(y_col = "observed",
      t_col = "arm", e_col = "e", m0_col = "m0", m1_col = "m1"),
    ratio = list(numerator_col = "y", denominator_col = "den"),
    clate = list(y_col = "observed", d_col = "arm", a_col = "arm", e_col = "e",
      m_y0_col = "m0", m_y1_col = "m1", m_d0_col = "d0", m_d1_col = "d1"))) else
    c(list(X = X), switch(family, np = list(Y = d$y), tdnn = list(Y = d$y),
      cate = list(pseudo_outcome = d$y), ratio = list(numerator = d$y, denominator = d$den),
      clate = list(Y = d$observed, D = arm, A = arm, e = d$e,
        m_y0 = d$m0, m_y1 = d$m1, m_d0 = d$d0, m_d1 = d$d1)))
  fn <- switch(family, np = np_dnn, tdnn = np_tdnn, cate = cate_dnn,
    ratio = affine_ratio_dnn, clate = clate_dnn)
  .with_localization_policy_warnings(do.call(fn, c(args, common)))
}

test_that("retained selectors distinguish missing fields from matrix NULL selectors", {
  targets <- rbind(c(.2, 5), c(10, .5))
  reordered <- data.frame(b = targets[, 2L], a = targets[, 1L])
  for (family in c("np", "tdnn", "cate", "ratio", "clate")) for (frame in c(FALSE, TRUE)) {
    fit <- .retained_field_fixture(family, frame)
    expect_true("x_cols" %in% names(fit$.model), info = paste(family, frame))
    if (!frame) expect_null(fit$.model[["x_cols", exact = TRUE]])
    else {
      expect_equal(predict(fit, newdata = reordered), predict(fit, X0 = targets))
      expect_gt(max(abs(predict(fit, X0 = targets) - predict(fit, X0 = targets[, 2:1]))), .1)
    }
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(predict(restored, X0 = targets, se.fit = TRUE), predict(fit, X0 = targets, se.fit = TRUE))
    expect_equal(diagnostics(restored), diagnostics(fit))
    bad <- fit
    bad$.model["x_cols"] <- NULL
    before <- serialize(bad, NULL)
    expect_error(predict(bad, newdata = reordered), "x_cols.*[Rr]efit")
    expect_error(diagnostics(bad), "x_cols.*[Rr]efit")
    expect_identical(serialize(bad, NULL), before)
    for (selectors in list(character(), NA_character_, "", 1:2, matrix(c("a", "b")))) {
      bad <- fit
      bad$.model$x_cols <- selectors
      expect_error(predict(bad, X0 = targets), "x_cols.*[Rr]efit")
    }
  }
})

test_that("retained selector validation preserves repeats and matrix-valued columns", {
  for (family in c("np", "tdnn", "cate", "ratio", "clate")) {
    for (selectors in list(c("a", "a"), "pair")) {
      fit <- .retained_field_fixture(family, frame = TRUE, x_cols = selectors)
      d <- data.frame(a = c(.2, 10), junk = c(900, 800))
      d$pair <- I(rbind(c(.2, 5), c(10, .5)))
      expected <- if (length(selectors) == 2L) cbind(d$a, d$a) else as.matrix(d$pair)
      expect_identical(fit$.model$x_cols, selectors)
      expect_equal(.with_localization_policy_warnings(predict(fit, newdata = d)),
        .with_localization_policy_warnings(predict(fit, X0 = expected)))
      expect_identical(.with_localization_policy_warnings(predict(unserialize(serialize(fit, NULL)), newdata = d, se.fit = TRUE)),
        .with_localization_policy_warnings(predict(fit, newdata = d, se.fit = TRUE)))
    }
  }
})

test_that("CATE diagnostics require current signal and overlap fields", {
  for (route in c("matrix", "frame_supplied", "frame_aipw")) {
    fit <- .retained_field_fixture("cate", frame = route != "matrix", supplied = route == "frame_supplied")
    for (field in c("psi", "t_vec", "e")) {
      bad <- fit
      bad$.model[field] <- NULL
      before <- serialize(bad, NULL)
      expect_error(diagnostics(bad), paste0(field, ".*[Rr]efit"))
      expect_error(predict(bad, X0 = fit$x0), paste0(field, ".*[Rr]efit"))
      expect_identical(serialize(bad, NULL), before)
    }
    if (route != "frame_aipw") {
      expect_true(all(c("t_vec", "e") %in% names(fit$.model)))
      expect_null(fit$.model$t_vec)
      expect_null(fit$.model$e)
      expect_false(any(diagnostics(fit)$overlap_available))
      expect_false(any(diagnostics(fit)$propensity_available))
    }
  }
  fit <- .retained_field_fixture("cate", frame = TRUE)
  for (field in c("psi", "t_vec", "e")) {
    for (value in list(NULL, 0, rep(NA_real_, 12L), matrix(fit$.model[[field]]))) {
      bad <- fit
      bad$.model[field] <- list(value)
      expect_error(diagnostics(bad), "[Rr]efit")
    }
  }
  for (field in c("t_vec", "e")) {
    bad <- fit
    bad$.model[[field]][1L] <- 2
    expect_error(diagnostics(bad), paste0(field, ".*[Rr]efit"))
  }
})

test_that("CLATE diagnostics retain required instrument record and training", {
  for (frame in c(FALSE, TRUE)) {
    fit <- .retained_field_fixture("clate", frame)
    expect_true(all(diagnostics(fit)$instrument_available))
    expect_identical(diagnostics(fit)$training, rep("supplied", 2L))
    paths <- list("iv", c("iv", "A"), c("iv", "training"), c("iv", "propensity"),
      c("iv", "propensity", "used"), c("iv", "propensity", "clipped"))
    for (path in paths) {
      bad <- fit
      if (length(path) == 1L) bad$.model[path] <- NULL
      if (length(path) == 2L) bad$.model[[path[1L]]][path[2L]] <- NULL
      if (length(path) == 3L) bad$.model[[path[1L]]][[path[2L]]][path[3L]] <- NULL
      before <- serialize(bad, NULL)
      expect_error(diagnostics(bad), "iv.*[Rr]efit")
      expect_error(predict(bad, X0 = fit$x0), "iv.*[Rr]efit")
      expect_identical(serialize(bad, NULL), before)
    }
    for (value in list(NULL, NA_character_, c("supplied", "oracle"), "retired")) {
      bad <- fit
      bad$.model$iv["training"] <- list(value)
      expect_error(diagnostics(bad), "training.*[Rr]efit")
    }
    for (field in c("A", "used", "clipped")) {
      bad <- fit
      if (field == "A") bad$.model$iv$A <- 1 else bad$.model$iv$propensity[[field]] <- 1
      expect_error(diagnostics(bad), "iv.*[Rr]efit")
    }
  }
})

test_that("overlap graphs reject lost retained state and preserve supplied NULL inputs", {
  skip_if_not_installed("ggplot2")
  fit <- .retained_field_fixture("cate", frame = TRUE)
  arm <- rep(c(0, 1), 6L)
  e <- seq(.15, .85, length.out = 12L)
  mass <- plot_overlap(fit, type = "mass", treatment = 1 - arm)
  expect_equal(mass$data$mass, c(as.numeric(crossprod(arm, fit$W)), as.numeric(crossprod(1 - arm, fit$W))))
  propensity <- plot_overlap(fit, type = "propensity", propensity = 1 - e)
  expect_equal(propensity$data$e, rep(e, ncol(fit$W)))
  for (field in c("psi", "t_vec", "e")) {
    bad <- fit
    bad$.model[field] <- NULL
    expect_error(plot_overlap(bad, type = "mass", treatment = 1 - arm), "[Rr]efit")
    expect_error(plot_overlap(bad, type = "propensity", propensity = 1 - e), "[Rr]efit")
  }
  supplied <- .retained_field_fixture("cate")
  expect_no_error(ggplot2::ggplot_build(plot_overlap(supplied, type = "mass", treatment = arm)))
  expect_no_error(.with_localization_policy_warnings(ggplot2::ggplot_build(plot_overlap(supplied, type = "propensity", propensity = e))))
  clate <- .retained_field_fixture("clate")
  clate$.model["iv"] <- NULL
  expect_error(plot_overlap(clate, type = "mass"), "iv.*[Rr]efit")
  expect_error(plot_overlap(clate, type = "propensity"), "iv.*[Rr]efit")
})
