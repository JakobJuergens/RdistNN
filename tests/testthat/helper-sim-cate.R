# Shared independent AIPW references for CATE simulation tests.
.test_aipw_score <- function(dat, e_hat, m0_hat, m1_hat) {
  pmax_e <- pmax(e_hat, 1e-8)
  pmax_1e <- pmax(1 - e_hat, 1e-8)
  m1_hat - m0_hat +
    dat$D * (dat$Y - m1_hat) / pmax_e -
    (1 - dat$D) * (dat$Y - m0_hat) / pmax_1e
}

.test_fit_aipw_nuisances <- function(dat, x_cols, train = rep(TRUE, nrow(dat)),
                                     clip = 1e-3) {
  e_fmla <- stats::as.formula(paste("D ~", paste(x_cols, collapse = " + ")))
  m_fmla <- stats::as.formula(paste("Y ~", paste(x_cols, collapse = " + ")))
  train_dat <- dat[train, , drop = FALSE]

  e_fit <- stats::glm(e_fmla, data = train_dat, family = stats::binomial())
  e_hat <- stats::predict(e_fit, newdata = dat, type = "response")
  e_hat <- pmin(pmax(e_hat, clip), 1 - clip)

  m0_fit <- stats::lm(m_fmla, data = train_dat[train_dat$D == 0L, , drop = FALSE])
  m1_fit <- stats::lm(m_fmla, data = train_dat[train_dat$D == 1L, , drop = FALSE])
  list(
    e_hat = e_hat,
    m0_hat = stats::predict(m0_fit, newdata = dat),
    m1_hat = stats::predict(m1_fit, newdata = dat)
  )
}
