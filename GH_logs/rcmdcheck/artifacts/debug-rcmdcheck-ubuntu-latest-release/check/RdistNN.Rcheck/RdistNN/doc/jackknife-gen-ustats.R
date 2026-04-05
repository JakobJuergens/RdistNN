## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(RdistNN)
set.seed(4)

## ----setup-example, eval=FALSE------------------------------------------------
# n <- 40
# p <- 2
# s <- 3L
# 
# X <- matrix(rnorm(n * p), n, p)
# Y <- rnorm(n)
# 
# # Toy compact-support example: build only the first s+1 ranked observations
# # for each target point. This illustrates the generic L-stat utilities below;
# # it is not the internal exact-DNN jackknife path used by np_dnn()/cate_dnn()/plqr_dnn().
# wt <- weights_dnn(
#   X = X, X0 = X,
#   method = "exact",
#   distance = "l2",
#   k = s, s = s,
#   return = "ord"
# )
# ranks <- t(wt$ord[seq_len(s + 1L), , drop = FALSE])  # T x (s+1)
# 
# # L-stat coefficients
# a <- rep(1, s)
# 
# # Scalar signal represented as n x 1 matrix
# phi <- cbind(Y)

## ----compare, eval=FALSE------------------------------------------------------
# jk_fast <- jackknife_fast_rank(ranks = ranks, a = a, phi = phi, normalize = TRUE)
# jk_slow <- jackknife_slow_rank(ranks = ranks, a = a, phi = phi, normalize = TRUE)
# 
# all.equal(jk_fast$theta0, jk_slow$theta0, tolerance = 1e-12)
# all.equal(jk_fast$theta_loo, jk_slow$theta_loo, tolerance = 1e-12)

## ----se, eval=FALSE-----------------------------------------------------------
# theta_loo <- jk_fast$theta_loo[, 1]
# theta_bar <- mean(theta_loo)
# n_ref <- length(theta_loo)
# 
# var_jk <- ((n_ref - 1) / n_ref) * sum((theta_loo - theta_bar)^2)
# se_jk <- sqrt(var_jk)
# se_jk

