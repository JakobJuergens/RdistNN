### Name: dnn_weights_exact_full_from_ranks
### Title: Exact DNN weights from ranks (full U-statistic, WOR)
### Aliases: dnn_weights_exact_full_from_ranks

### ** Examples

set.seed(1)
n <- 10
s <- 4
k <- 2
ranks <- sample.int(n) # ranks[i] in 1..n
out <- dnn_weights_exact_full_from_ranks(ranks, k = k, s = s)
sum(out$w) # 1
sum(out$w_rank) # 1



