# Visual regression baselines for the plotting panels.
# These snapshots are opt-in: SVG rendering differs across
# platforms and font stacks, so they run only when RDISTNN_VDIFFR=true is
# set explicitly (a deliberate local run), never on CI or CRAN.
#
# Regenerate baselines with:
#   RDISTNN_VDIFFR=true Rscript -e 'testthat::test_local(filter = "vdiffr")'
# then review and commit tests/testthat/_snaps/plot-vdiffr/.

.vdiffr_active <- function() {
  identical(Sys.getenv("RDISTNN_VDIFFR"), "true") &&
    requireNamespace("vdiffr", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)
}

test_that("visual baselines for the core panels", {
  skip_on_cran()
  skip_if(!.vdiffr_active(),
          "visual baselines run only with RDISTNN_VDIFFR=true")

  set.seed(2026)
  n <- 120L
  df <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  df$e <- plogis(0.5 * df$X1)
  df$T <- as.integer(runif(n) < df$e)
  df$m0 <- df$X1
  df$m1 <- df$X1 + 1
  df$Y <- ifelse(df$T == 1L, df$m1, df$m0) + rnorm(n, sd = 0.3)
  X0 <- cbind(seq(-1, 1, length.out = 6L), 0)

  fit <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    e_col = "e", m0_col = "m0", m1_col = "m1",
    X0 = X0, k = 1L, s = 15L, se = TRUE, return = c("W", "dx")
  )

  vdiffr::expect_doppelganger("cate estimates profile", plot_estimates(fit))
  vdiffr::expect_doppelganger("weights kernel", plot_weights(fit, type = "kernel", targets = 3L))
  vdiffr::expect_doppelganger("weights cumulative", plot_weights(fit, type = "cumulative"))
  vdiffr::expect_doppelganger("overlap mass", plot_overlap(fit, type = "mass"))
  vdiffr::expect_doppelganger("overlap propensity", plot_overlap(fit, type = "propensity", targets = 3L))
  vdiffr::expect_doppelganger("diagnostics facets", ggplot2::autoplot(diagnostics(fit)))
  vdiffr::expect_doppelganger("influence", plot_influence(fit, targets = 3L))

  grid_tuning <- new_tuning(
    name = "vdiffr_grid",
    eval = function(...) list(s = 15L, meta = list(lambda = 1)),
    eval_grid = function(...) data.frame(
      lambda = c(0.75, 1, 1.25),
      s = c(10L, 15L, 20L),
      stringsAsFactors = FALSE
    )
  )
  fit_grid <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    e_col = "e", m0_col = "m0", m1_col = "m1",
    X0 = X0, k = 1L, tuning = grid_tuning, tuning_mode = "grid",
    se = TRUE, return = c("W", "dx", "ranks")
  )
  vdiffr::expect_doppelganger("grid candidate estimates", plot_estimates(fit_grid, candidates = "all"))
  vdiffr::expect_doppelganger("grid candidate diagnostics", ggplot2::autoplot(diagnostics(fit_grid, candidates = "all")))
  vdiffr::expect_doppelganger("grid tuning sweep", plot_tuning(fit_grid))

  demo <- plqr_demo_data(n = 90L, p = 2L, p_d = 2L, tau = 0.5, seed = 7)
  dfp <- demo$data
  fit_plqr <- plqr_dnn(
    X = as.matrix(dfp[paste0("X", 1:2)]),
    D = as.matrix(dfp[paste0("D", 1:2)]),
    g = as.numeric(dfp$g_oracle),
    mD = as.matrix(dfp[paste0("mD", 1:2)]),
    Y = as.numeric(dfp$Y),
    theta_eval = demo$truth$theta_true, tau = 0.5,
    X0 = cbind(seq(-0.5, 0.5, length.out = 5L), 0),
    k = 1L, s = 10L, se = TRUE
  )
  vdiffr::expect_doppelganger("plqr moments", ggplot2::autoplot(fit_plqr))
  vdiffr::expect_doppelganger("plqr conditioning", plot_conditioning(fit_plqr))
  vdiffr::expect_doppelganger("plqr calibration", plot_calibration(fit_plqr))

  set.seed(2027)
  x_mat <- matrix(rnorm(100L * 2L), 100L, 2L,
                  dimnames = list(NULL, c("X1", "X2")))
  y_vec <- x_mat[, 1L] + rnorm(100L, sd = 0.3)
  fit_tdnn <- np_tdnn(
    X = x_mat, X0 = cbind(seq(-1, 1, length.out = 6L), 0), Y = y_vec,
    s1 = 5L, s2 = 10L, se = FALSE
  )
  vdiffr::expect_doppelganger("tdnn scales", plot_tdnn_scales(fit_tdnn))
})
