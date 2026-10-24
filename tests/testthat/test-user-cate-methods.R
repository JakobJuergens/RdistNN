.make_cate_fit <- function(n = 40L, p = 2L, m = 3L, se = TRUE, seed = 456L) {
  .test_rng_scope({
    X <- matrix(rnorm(n * p), n, p)
    X0 <- matrix(rnorm(m * p), m, p)
    psi <- X[, 1] - 0.25 * X[, min(2L, p)] + rnorm(n, sd = 0.2)
    fit <- cate_dnn(
      X = X, X0 = X0, pseudo_outcome = psi,
      distance = "l2", method = "exact",
      k = 2L, s = 3L, se = se, return = "W"
    )
    list(fit = fit, X0 = X0)
  }, seed = seed)
}

test_that("summary.rdistnn_cate_result returns expected structure", {
  obj <- .make_cate_fit(se = TRUE)
  s <- summary(obj$fit)

  expect_s3_class(s, "rdistnn_summary_cate_result")
  expect_true(is.data.frame(s$coefficients))
  expect_equal(
    colnames(s$coefficients),
    c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  )
  expect_equal(nrow(s$coefficients), length(obj$fit$tau))
  expect_true(all(is.finite(s$coefficients$estimate)))
})

test_that("summary handles missing se by returning NA inferential columns", {
  obj <- .make_cate_fit(se = FALSE)
  s <- summary(obj$fit)

  expect_false(isTRUE(s$has_se))
  expect_true(all(is.na(s$coefficients$std.error)))
  expect_true(all(is.na(s$coefficients$statistic)))
  expect_true(all(is.na(s$coefficients$p.value)))
})

test_that("confint.rdistnn_cate_result returns valid intervals and honors parm", {
  obj <- .make_cate_fit(se = TRUE)
  ci_all <- confint(obj$fit)
  ci_one <- confint(obj$fit, parm = 2L)

  expect_true(is.matrix(ci_all))
  expect_equal(dim(ci_all), c(length(obj$fit$tau), 2L))
  expect_true(all(ci_all[, "conf.low"] <= ci_all[, "conf.high"]))
  expect_equal(dim(ci_one), c(1L, 2L))
  expect_true(grepl("target_2", rownames(ci_one)[1], fixed = TRUE))
})

test_that("confint errors when standard errors are not available", {
  obj <- .make_cate_fit(se = FALSE)
  missing_extra <- expect_error(
    confint(obj$fit),
    "se = TRUE",
    ignore.case = TRUE, class = "rdistnn_extra_unavailable"
  )
  expect_identical(missing_extra$code, "extra_unavailable")
  expect_match(conditionMessage(missing_extra), "[Rr]efit.*se = TRUE")
})

test_that("predict.rdistnn_cate_result reproduces fitted tau on same X0", {
  obj <- .make_cate_fit(se = TRUE)
  pred <- predict(obj$fit, X0 = obj$X0)
  expect_equal(pred, obj$fit$tau, tolerance = 1e-12)
})

test_that("predict supports se.fit and confidence intervals", {
  obj <- .make_cate_fit(se = TRUE)

  pred_se <- predict(obj$fit, X0 = obj$X0, se.fit = TRUE)
  expect_true(is.list(pred_se))
  expect_true(all(c("fit", "se.fit") %in% names(pred_se)))
  expect_equal(length(pred_se$fit), nrow(obj$X0))
  expect_equal(length(pred_se$se.fit), nrow(obj$X0))
  expect_true(all(pred_se$se.fit >= 0))

  pred_ci <- predict(obj$fit, X0 = obj$X0, interval = "confidence")
  expect_true(is.matrix(pred_ci))
  expect_equal(dim(pred_ci), c(nrow(obj$X0), 3L))
  expect_equal(colnames(pred_ci), c("fit", "lwr", "upr"))
  expect_true(all(pred_ci[, "lwr"] <= pred_ci[, "fit"]))
  expect_true(all(pred_ci[, "fit"] <= pred_ci[, "upr"]))

  pred_ci_se <- predict(obj$fit, X0 = obj$X0, interval = "confidence", se.fit = TRUE)
  expect_true(is.list(pred_ci_se))
  expect_true(all(c("fit", "se.fit") %in% names(pred_ci_se)))
  expect_true(is.matrix(pred_ci_se$fit))
  expect_equal(dim(pred_ci_se$fit), c(nrow(obj$X0), 3L))
  expect_equal(length(pred_ci_se$se.fit), nrow(obj$X0))
})

test_that("predict accepts single x0 and validates target inputs", {
  obj <- .make_cate_fit(se = TRUE)
  x0 <- obj$X0[1, ]
  pred <- predict(obj$fit, x0 = x0)
  expect_true(is.numeric(pred))
  expect_length(pred, 1L)

  expect_error(
    predict(obj$fit),
    "Provide `newdata`, `X0`, or a single `x0`",
    ignore.case = TRUE
  )
  expect_error(
    predict(obj$fit, X0 = obj$X0, x0 = x0),
    "either `X0` or `x0`",
    ignore.case = TRUE
  )
  expect_error(
    predict(obj$fit, newdata = obj$X0, X0 = obj$X0),
    "either `newdata` or `X0`/`x0`",
    fixed = TRUE
  )
})

test_that("predict handles data.frame newdata via stored x_cols", {
  d <- .make_df_data(n = 55L, p = 2L, m = 4L, seed = 31L)
  fit <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D", m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = TRUE
  )

  nd <- data.frame(
    X1 = d$X0[, 1],
    X2 = d$X0[, 2],
    extra = seq_len(nrow(d$X0))
  )
  pred_df <- predict(fit, newdata = nd)
  pred_X0 <- predict(fit, X0 = as.matrix(nd[, d$x_cols, drop = FALSE]))

  expect_equal(pred_df, pred_X0, tolerance = 1e-12)
})

test_that("predict validates se.fit and target values before weight computation", {
  obj <- .make_cate_fit(se = TRUE)

  expect_error(
    predict(obj$fit, X0 = obj$X0, se.fit = "yes"),
    "se.fit.*TRUE or FALSE"
  )
  expect_error(
    predict(obj$fit, X0 = obj$X0, se.fit = NA),
    "se.fit.*TRUE or FALSE"
  )
  missing_target <- expect_error(
    predict(obj$fit, x0 = c(NA_real_, 0)),
    "x0 cannot contain `NA` values",
    fixed = TRUE, class = "rdistnn_input_error"
  )
  expect_identical(missing_target$code, "argument_invalid")
  expect_identical(missing_target$argument, "x0")

  bad_X0 <- obj$X0
  bad_X0[1L, 1L] <- Inf
  nonfinite_target <- expect_error(
    predict(obj$fit, X0 = bad_X0),
    "X0 cannot contain non-finite values",
    fixed = TRUE, class = "rdistnn_input_error"
  )
  expect_identical(nonfinite_target$code, "argument_invalid")
  expect_identical(nonfinite_target$argument, "X0")

  bad_newdata <- obj$X0
  bad_newdata[1L, 1L] <- NaN
  expect_error(
    predict(obj$fit, newdata = bad_newdata),
    "newdata cannot contain non-finite values",
    fixed = TRUE
  )
  expect_error(
    predict(obj$fit, newdata = data.frame(X1 = "a", X2 = "b")),
    "newdata must be numeric",
    fixed = TRUE
  )
  expect_error(
    predict(obj$fit, X0 = matrix(numeric(0), nrow = 0L, ncol = ncol(obj$X0))),
    "X0 must have at least 1 row",
    fixed = TRUE
  )
})

test_that("predict validates data-frame newdata columns", {
  d <- .make_df_data(n = 55L, p = 2L, m = 4L, seed = 34L)
  fit <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D", m0_col = "m0_oracle", m1_col = "m1_oracle", e_col = "e_oracle",
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_error(
    predict(fit, newdata = data.frame(X1 = d$X0[, 1L])),
    "X2.*x_cols"
  )

  nd <- data.frame(X1 = d$X0[, 1L], X2 = d$X0[, 2L])
  nd$X2[1L] <- NA_real_
  expect_error(
    predict(fit, newdata = nd),
    "newdata cannot contain `NA` values",
    fixed = TRUE
  )
})

test_that("predict validates fitted model state before computing", {
  obj <- .make_cate_fit(se = FALSE)

  bad_psi <- obj$fit
  bad_psi$.model$psi <- bad_psi$.model$psi[-1L]
  expect_error(
    predict(bad_psi, X0 = obj$X0),
    "object\\$\\.model\\$psi.*[Rr]efit"
  )

  bad_x <- obj$fit
  bad_x$.model$x_mat[1L, 1L] <- Inf
  error <- expect_error(predict(bad_x, X0 = obj$X0), class = "rdistnn_state_error")
  expect_identical(error$code, "state_invalid")
  expect_match(conditionMessage(error), "[Rr]efit")

  bad_k <- obj$fit
  bad_k$.model$k <- bad_k$.model$s + 1L
  error <- expect_error(
    predict(bad_k, X0 = obj$X0),
    "meta\\$k/\\.model\\$k.*inconsistent.*[Rr]efit",
    class = "rdistnn_state_error"
  )
  expect_identical(error$code, "state_invalid")

  bad_s <- obj$fit
  bad_s$.model$s <- nrow(bad_s$.model$x_mat) + 1L
  expect_error(
    predict(bad_s, X0 = obj$X0),
    "meta\\$s/\\.model\\$s.*inconsistent.*[Rr]efit",
    class = "rdistnn_state_error"
  )

  bad_method <- obj$fit
  bad_method$.model$method <- "bogus"
  expect_error(
    predict(bad_method, X0 = obj$X0),
    "method.*exact.*asymp.*[Rr]efit",
    class = "rdistnn_state_error"
  )

  bad_normalize <- obj$fit
  bad_normalize$.model$normalize <- NA
  expect_error(
    predict(bad_normalize, X0 = obj$X0),
    "normalize.*TRUE or FALSE.*[Rr]efit",
    class = "rdistnn_state_error"
  )
})

test_that("predict interval requests enforce jackknife feasibility", {
  obj <- .make_cate_fit(se = FALSE)
  model <- obj$fit$.model
  full_order <- cate_dnn(X = model$x_mat, X0 = obj$X0, pseudo_outcome = model$psi,
    k = model$k, s = nrow(model$x_mat), se = FALSE, return = character())

  expect_error(
    predict(full_order, X0 = obj$X0, se.fit = TRUE),
    "s < n",
    fixed = TRUE
  )
  expect_error(
    predict(full_order, X0 = obj$X0, interval = "confidence"),
    "s < n",
    fixed = TRUE
  )
})

test_that("methods validate level/parm and error on missing model state", {
  obj <- .make_cate_fit(se = TRUE)

  expect_error(summary(obj$fit, level = 0), "numeric scalar in \\(0,1\\)", ignore.case = TRUE)
  expect_error(confint(obj$fit, level = 1), "numeric scalar in \\(0,1\\)", ignore.case = TRUE)
  expect_error(
    confint(obj$fit, parm = c(0L, 2L)),
    "valid target indices",
    ignore.case = TRUE
  )
  expect_error(
    confint(obj$fit, parm = 1.5),
    "integer target indices",
    ignore.case = TRUE
  )
  expect_error(
    confint(obj$fit, parm = Inf),
    "non-finite",
    ignore.case = TRUE
  )
  expect_error(
    confint(obj$fit, parm = 1e20),
    "valid target indices",
    ignore.case = TRUE
  )
  expect_error(
    confint(obj$fit, parm = "target_1"),
    "valid target indices",
    ignore.case = TRUE
  )
  expect_error(
    predict(obj$fit, X0 = obj$X0, interval = "confidence", level = 1.5),
    "numeric scalar in \\(0,1\\)",
    ignore.case = TRUE
  )

  fit_no_model <- obj$fit
  fit_no_model$.model <- NULL
  error <- expect_error(predict(fit_no_model, X0 = obj$X0), class = "rdistnn_state_error")
  expect_identical(error$code, "state_invalid")
  expect_match(error$argument, "object$.model", fixed = TRUE)
  expect_match(conditionMessage(error), "[Rr]efit")
})

test_that("print.rdistnn_summary_cate_result runs", {
  obj <- .make_cate_fit(se = TRUE)
  s <- summary(obj$fit)
  txt_summary <- capture.output(print(s))
  expect_true(any(grepl("<rdistnn_summary_cate_result>", txt_summary, fixed = TRUE)))
})

test_that("print.rdistnn_cate_result reports standard-error availability invisibly", {
  with_se <- .make_cate_fit(se = TRUE)$fit
  with_se_text <- capture.output(with_se_returned <- print(with_se))

  expect_identical(with_se_returned, with_se)
  expect_true(any(grepl("<rdistnn_cate_result>", with_se_text, fixed = TRUE)))
  expect_true(any(grepl("SEs available", with_se_text, fixed = TRUE)))

  without_se <- .make_cate_fit(se = FALSE)$fit
  without_se_text <- capture.output(without_se_returned <- print(without_se))

  expect_identical(without_se_returned, without_se)
  expect_true(any(grepl("SEs not requested", without_se_text, fixed = TRUE)))
})

test_that("coef.rdistnn_cate_result returns named target estimates", {
  obj <- .make_cate_fit(se = FALSE)
  expect_equal(unname(coef(obj$fit)), obj$fit$tau)
  expect_identical(names(coef(obj$fit)), paste0("target_", seq_along(obj$fit$tau)))
})

test_that("two-target CATE coefficients retain target identities", {
  fit <- cate_dnn(X = matrix(1:8), pseudo_outcome = 1:8, X0 = matrix(c(0, 9)), s = 3, se = FALSE)
  expect_identical(names(coef(fit)), c("target_1", "target_2"))
})
