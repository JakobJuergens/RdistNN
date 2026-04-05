skip_if_not(
  all(vapply(
    c("cate_dnn", "summary.cate_result", "confint.cate_result", "predict.cate_result"),
    exists, logical(1),
    mode = "function"
  )),
  "cate_result methods or required helpers not found; skipping file"
)

set.seed(456)

.make_cate_fit <- function(n = 40L, p = 2L, m = 3L, se = TRUE) {
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  psi <- X[, 1] - 0.25 * X[, min(2L, p)] + rnorm(n, sd = 0.2)
  fit <- cate_dnn(
    X = X, X0 = X0, pseudo_outcome = psi,
    distance = "l2", method = "exact",
    k = 2L, s = 3L, se = se, return = "W"
  )
  list(fit = fit, X0 = X0)
}

test_that("summary.cate_result returns expected structure", {
  obj <- .make_cate_fit(se = TRUE)
  s <- summary(obj$fit)

  expect_s3_class(s, "summary.cate_result")
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

test_that("confint.cate_result returns valid intervals and honors parm", {
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
  expect_error(
    confint(obj$fit),
    "se = TRUE",
    ignore.case = TRUE
  )
})

test_that("predict.cate_result reproduces fitted tau on same X0", {
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
})

test_that("predict handles data.frame newdata via stored x_cols", {
  d <- .make_df_data(n = 55L, p = 2L, m = 4L, seed = 31L)
  fit <- cate_dnn(
    data = d$df, x_cols = d$x_cols,
    y_col = "Y", t_col = "D", m_col = "m_oracle", e_col = "e_oracle",
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
    predict(obj$fit, X0 = obj$X0, interval = "confidence", level = 1.5),
    "numeric scalar in \\(0,1\\)",
    ignore.case = TRUE
  )

  fit_no_model <- obj$fit
  fit_no_model$.model <- NULL
  expect_error(
    predict(fit_no_model, X0 = obj$X0),
    "does not contain fitted model state",
    ignore.case = TRUE
  )
})

test_that("print.summary.cate_result and legacy print.cate_dnn_result run", {
  obj <- .make_cate_fit(se = TRUE)
  s <- summary(obj$fit)
  txt_summary <- capture.output(print(s))
  expect_true(any(grepl("<summary.cate_result>", txt_summary, fixed = TRUE)))

  legacy <- list(
    tau = obj$fit$tau,
    se = obj$fit$se,
    meta = list(
      n = obj$fit$meta$n,
      p = obj$fit$meta$p,
      m = obj$fit$meta$m,
      k = obj$fit$meta$k,
      s = obj$fit$meta$s,
      method = obj$fit$meta$method,
      weight_name = "dnn_exact",
      distance_name = obj$fit$meta$distance
    )
  )
  class(legacy) <- "cate_dnn_result"
  txt_legacy <- capture.output(print(legacy))
  expect_true(any(grepl("<cate_dnn_result>", txt_legacy, fixed = TRUE)))
})
