## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(RdistNN)
set.seed(1)

## ----np, eval=FALSE-----------------------------------------------------------
# demo_np <- np_demo_data(n = 250, p = 2, design = "sinusoidal", seed = 1)
# 
# fit_np <- np_dnn(
#   data = demo_np$data,
#   x_cols = c("X1", "X2"),
#   y_col = "Y",
#   x0 = c(0, 0),
#   method = "exact",
#   k = 2L, s = 3L,
#   se = TRUE,
#   jackknife = "fast",
#   return = c("W", "ess")
# )
# 
# fit_np$mu
# fit_np$se

## ----cate, eval=FALSE---------------------------------------------------------
# demo_cate <- cate_demo_data(n = 300, p = 2, design = "hetero", seed = 1)
# df <- demo_cate$data
# df$m1_oracle <- df$m0_oracle + df$tau_oracle
# 
# fit_cate <- cate_dnn(
#   data = df,
#   x_cols = c("X1", "X2"),
#   y_col = "Y",
#   t_col = "D",
#   m0_col = "m0_oracle",
#   m1_col = "m1_oracle",
#   e_col = "e_oracle",
#   x0 = c(0, 0),
#   method = "exact",
#   k = 2L, s = 3L,
#   se = TRUE,
#   return = c("psi", "nuisance")
# )
# 
# summary(fit_cate)
# confint(fit_cate)
# predict(fit_cate, x0 = c(0.2, -0.1), interval = "confidence")

## ----plqr, eval=FALSE---------------------------------------------------------
# demo_plqr <- plqr_demo_data(n = 250, p = 2, p_d = 2, tau = 0.5, seed = 1)
# 
# fit_plqr <- plqr_dnn(
#   data = demo_plqr$data,
#   x_cols = c("X1", "X2"),
#   D_cols = c("D1", "D2"),
#   q_col = "q_oracle",
#   mD_cols = c("mD1", "mD2"),
#   y_col = "Y",
#   beta = demo_plqr$truth$beta_true,
#   tau = 0.5,
#   x0 = c(0, 0),
#   method = "exact",
#   k = 2L, s = 3L,
#   se = TRUE,
#   return = c("g", "ess")
# )
# 
# fit_plqr$moment
# fit_plqr$se

