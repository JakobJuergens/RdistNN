## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(RdistNN)
set.seed(2)

## ----basic, eval=FALSE--------------------------------------------------------
# n <- 150
# p <- 2
# X <- matrix(rnorm(n * p), n, p)
# X0 <- matrix(rnorm(8), 4, p)
# 
# wt_exact <- weights_dnn(
#   X = X, X0 = X0,
#   distance = "l2",
#   method = "exact",
#   k = 2L, s = 4L,
#   return = c("W", "ranks", "ord", "dx")
# )
# 
# wt_asymp <- weights_dnn(
#   X = X, X0 = X0,
#   distance = "l2",
#   method = "asymp",
#   k = 2L, s = 4L,
#   normalize = TRUE,
#   return = "W"
# )
# 
# colSums(wt_exact$W)
# colSums(wt_asymp$W)

## ----single-target, eval=FALSE------------------------------------------------
# x0 <- c(0, 0)
# wt_single <- weights_dnn(
#   X = X, x0 = x0,
#   method = "exact",
#   k = 2L, s = 4L,
#   return = c("W", "ranks")
# )
# 
# dim(wt_single$W)  # n x 1

## ----compare, eval=FALSE------------------------------------------------------
# mu_signal <- rnorm(n)
# mu_exact <- as.numeric(crossprod(mu_signal, wt_exact$W))
# mu_asymp <- as.numeric(crossprod(mu_signal, wt_asymp$W))
# 
# cbind(mu_exact, mu_asymp)

