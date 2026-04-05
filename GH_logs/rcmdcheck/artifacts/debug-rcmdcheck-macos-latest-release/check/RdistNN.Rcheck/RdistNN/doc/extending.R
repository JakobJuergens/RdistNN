## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(RdistNN)
set.seed(3)

## ----distance-simple, eval=FALSE----------------------------------------------
# # L-infinity distance
# linf <- as_distance(
#   function(X, x0, ctx = NULL) {
#     apply(abs(X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)), 1L, max)
#   },
#   name = "linf"
# )
# 
# register_distance("linf", linf)
# distance_names()

## ----distance-eval, eval=FALSE------------------------------------------------
# X <- matrix(rnorm(40), 20, 2)
# x0 <- c(0, 0)
# 
# distance_eval(linf, X, x0)

## ----distance-batch, eval=FALSE-----------------------------------------------
# scaled_l2 <- new_distance(
#   name = "scaled_l2_demo",
#   setup = function(X) {
#     scale <- apply(X, 2L, stats::sd)
#     scale[!is.finite(scale) | scale <= 0] <- 1
#     list(scale = scale)
#   },
#   eval_batch = function(X, X0, ctx) {
#     Xs <- sweep(X, 2L, ctx$scale, "/")
#     X0s <- sweep(X0, 2L, ctx$scale, "/")
# 
#     vapply(seq_len(nrow(X0s)), function(j) {
#       rowSums((Xs - matrix(X0s[j, ], nrow(Xs), ncol(Xs), byrow = TRUE))^2)
#     }, numeric(nrow(Xs)))
#   }
# )
# 
# X <- matrix(rnorm(80), 40, 2)
# X0 <- matrix(rnorm(12), 6, 2)
# 
# ctx_dist <- scaled_l2$setup(X)
# D <- distance_eval_batch(scaled_l2, X, X0, ctx = ctx_dist)
# dim(D)

## ----rank-simple, eval=FALSE--------------------------------------------------
# rw_inv <- as_rank_weight(
#   function(ranks, ctx = NULL) {
#     w <- 1 / as.numeric(ranks)
#     w / sum(w)
#   },
#   name = "inv_rank"
# )
# 
# register_rank("inv_rank", rw_inv)
# rank_names()

## ----rank-support, eval=FALSE-------------------------------------------------
# make_top_m <- function(m) {
#   stopifnot(length(m) == 1L, m >= 1L)
# 
#   as_rank_weight(
#     function(ranks) {
#       keep <- as.numeric(ranks <= m)
#       keep / sum(keep)
#     },
#     name = paste0("top_", m)
#   )
# }
# 
# rw_top10 <- make_top_m(10L)

## ----rank-batch, eval=FALSE---------------------------------------------------
# rw_exp <- new_rank_weight(
#   name = "exp_decay",
#   eval_batch = function(R, ctx = NULL) {
#     lambda <- if (is.null(ctx) || is.null(ctx$lambda)) 0.25 else ctx$lambda
#     W <- exp(-lambda * (R - 1L))
#     sweep(W, 2L, colSums(W), "/")
#   }
# )
# 
# R <- cbind(1:5, c(2L, 1L, 3L, 4L, 5L))
# rank_weight_eval_batch(rw_exp, R, ctx = list(lambda = 0.5))

## ----weights-matrix, eval=FALSE-----------------------------------------------
# n <- 80
# p <- 2
# X <- matrix(rnorm(n * p), n, p)
# X0 <- matrix(rnorm(4 * p), 4, p)
# 
# # If only the distance needs setup, leave ctx = NULL and let the distance
# # auto-setup path run.
# out1 <- weights_matrix(
#   X = X,
#   X0 = X0,
#   ranking = scaled_l2,
#   weight = rw_top10,
#   return = c("W", "dx", "ord", "ranks")
# )
# 
# # If both the distance and the weight need context, build one shared list.
# ctx_both <- c(scaled_l2$setup(X), list(lambda = 0.35))
# out2 <- weights_matrix(
#   X = X,
#   X0 = X0,
#   ranking = scaled_l2,
#   weight = rw_exp,
#   ctx = ctx_both,
#   return = c("W", "dx", "ord", "ranks")
# )
# 
# dim(out2$W)

## ----weights-debug, eval=FALSE------------------------------------------------
# out_dbg <- weights_matrix(
#   X = X,
#   X0 = X0,
#   ranking = scaled_l2,
#   weight = rw_exp,
#   ctx = ctx_both,
#   return = c("dx", "ord", "ranks", "W")
# )
# 
# head(out_dbg$dx[, 1])
# head(out_dbg$ord[, 1])
# head(out_dbg$ranks[, 1])
# head(out_dbg$W[, 1])

## ----score-simple, eval=FALSE-------------------------------------------------
# sc_avg2 <- as_score(
#   function(idxs, data, ctx = NULL, rng = NULL) {
#     mean(data$Y[idxs])
#   },
#   s = 2L,
#   name = "avg2"
# )
# 
# register_score(sc_avg2)
# list_scores()
# 
# df <- data.frame(Y = rnorm(10))
# subs <- rbind(c(1L, 2L), c(4L, 8L), c(3L, 9L))
# score_eval_batch(sc_avg2, subs = subs, data = df)

## ----score-setup, eval=FALSE--------------------------------------------------
# sc_pair_gap <- new_score(
#   name = "pair_gap_sq",
#   order = 2L,
#   setup = function(data, ...) list(y = as.numeric(data$Y)),
#   eval_batch = function(subs, data, ctx, rng = NULL) {
#     (ctx$y[subs[, 1L]] - ctx$y[subs[, 2L]])^2
#   }
# )
# 
# df <- data.frame(Y = rnorm(10))
# ctx_score <- score_setup(sc_pair_gap, df)
# subs <- rbind(c(1L, 2L), c(2L, 5L), c(3L, 4L))
# 
# score_eval_batch(sc_pair_gap, subs, data = df, ctx = ctx_score)

## ----score-diff, eval=FALSE---------------------------------------------------
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

## ----estimator-pattern, eval=FALSE--------------------------------------------
# my_local_mean <- function(X, x0, y, ranking = "l2", weight) {
#   W <- weights_matrix(
#     X = X,
#     x0 = x0,
#     ranking = ranking,
#     weight = weight,
#     return = "W"
#   )$W[, 1L]
# 
#   sum(W * y)
# }

