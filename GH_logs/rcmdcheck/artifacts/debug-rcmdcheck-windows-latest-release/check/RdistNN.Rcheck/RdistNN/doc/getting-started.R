## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(RdistNN)

## ----weights, eval=FALSE------------------------------------------------------
# set.seed(1)
# n <- 120
# p <- 2
# X <- matrix(rnorm(n * p), n, p)
# X0 <- matrix(rnorm(6), 3, p)
# 
# # Exact finite-sample DNN weights
# wt_exact <- weights_dnn(
#   X = X, X0 = X0,
#   distance = "l2",
#   method = "exact",
#   k = 2L, s = 3L,
#   return = c("W", "ranks", "ord")
# )
# 
# # Asymptotic DNN weights
# wt_asymp <- weights_dnn(
#   X = X, X0 = X0,
#   distance = "l2",
#   method = "asymp",
#   k = 2L, s = 3L,
#   normalize = TRUE,
#   return = "W"
# )
# 
# dim(wt_exact$W)

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
#   return = c("W", "ess", "neighbors")
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

## ----custom, eval=FALSE-------------------------------------------------------
# register_distance("linf_demo", function(X, x0, ctx = NULL) {
#   apply(abs(X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)), 1L, max)
# })
# 
# register_rank("inv_rank_demo", function(ranks) {
#   w <- 1 / as.numeric(ranks)
#   w / sum(w)
# })
# 
# X <- as.matrix(np_demo_data(n = 80, p = 2, seed = 1)$data[, c("X1", "X2")])
# 
# out <- weights_matrix(
#   X = X,
#   x0 = c(0, 0),
#   ranking = "linf_demo",
#   weight = get_rank("inv_rank_demo"),
#   return = c("W", "ranks")
# )
# dim(out$W)

## ----custom-score, eval=FALSE-------------------------------------------------
# sc_affine <- new_score(
#   name = "affine_order1",
#   order = 1L,
#   eval_param = function(idxs, data, theta, ctx, rng = NULL) {
#     theta[1] * data$Y[idxs] + theta[2]
#   },
#   grad = function(idxs, data, theta, ctx, rng = NULL) {
#     c(data$Y[idxs], 1)
#   }
# )
# 
# df <- data.frame(Y = c(1, 2, 4))
# path <- score_diff_path(sc_affine, data = df, method = "auto")
# 
# path$value(2L, theta = c(0.5, 1))
# path$grad(2L, theta = c(0.5, 1))

