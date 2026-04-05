### Name: dnn_weights_asymp_full_from_ranks
### Title: Asymptotic DNN weights from ranks (full U-statistic, WOR)
### Aliases: dnn_weights_asymp_full_from_ranks

### ** Examples

set.seed(1)
n <- 50
s <- 10
k <- 3
ranks <- sample.int(n)
out <- dnn_weights_asymp_full_from_ranks(ranks, k = k, s = s)
sum(out$w) # 1
sum(out$w_rank) # 1



