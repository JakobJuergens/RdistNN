test_that("engine + CATE score agree (compiled vs R fallback)", {
  skip_if_not(exists("u_statistic", mode = "function"))
  skip_if_not(has_cate_cpp(), "C++ CATE functions not available")

  demo <- cate_demo_data(n = 180, p = 2, seed = 77)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc <- score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho = rho, localizer = "dnnk", k = 5)

  # Compiled path
  u_cpp <- u_statistic(data = df, score = sc, sampling = list(mode = "fixedN", N = 25000))
  fit_cpp <- fit(u_cpp)$value

  # Force R fallback by nulling ptr in the stored ctx *inside* the fit object if available,
  # or re-run with a score that has no compiled ctx.
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL
  # If your u_statistic() accepts a prebuilt ctx, use it here; otherwise emulate via a score wrapper:
  sc_r <- sc
  sc_r$setup <- function(data, ...) { ctx }  # ensures no ptr
  u_r <- u_statistic(data = df, score = sc_r, sampling = list(mode = "fixedN", N = 25000))
  fit_r <- fit(u_r)$value

  expect_equal(fit_cpp, fit_r, tolerance = 1e-8)
})
