pkgname <- "RdistNN"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "RdistNN-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('RdistNN')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("as_distance")
### * as_distance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as_distance
### Title: Wrap a function as a distance spec
### Aliases: as_distance

### ** Examples

spec <- as_distance(
  function(X, x0, ctx = NULL) {
    rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  }
)
x_mat <- matrix(1:9, nrow = 3)
distance_eval(spec, x_mat, c(0, 0, 0))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as_distance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("as_rank_weight")
### * as_rank_weight

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as_rank_weight
### Title: Wrap a function as a rank_weight spec
### Aliases: as_rank_weight

### ** Examples

rw <- as_rank_weight(function(ranks) {
  w <- 1 / as.numeric(ranks)
  w / sum(w)
})
rank_weight_eval(rw, c(3L, 1L, 2L))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as_rank_weight", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("as_score")
### * as_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as_score
### Title: Coerce a plain function to a score (simple cases)
### Aliases: as_score

### ** Examples

sc <- as_score(function(idxs, data) prod(data[idxs]), s = 2L, name = "prod2")
score_eval(sc, c(1L, 3L), c(2, 4, 6), NULL)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cate_demo_data")
### * cate_demo_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cate_demo_data
### Title: Minimal CATE demo data (oracle nuisances + pseudo-outcome)
### Aliases: cate_demo_data

### ** Examples

demo <- cate_demo_data(n = 200, p = 2, design = "hetero", seed = 1)
str(demo$data)
demo$tau_x0




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cate_demo_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cate_demo_subsets")
### * cate_demo_subsets

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cate_demo_subsets
### Title: Quick helper: sample a batch of subsets (for U-stat experiments)
### Aliases: cate_demo_subsets

### ** Examples

subs <- cate_demo_subsets(n = 200, K = 1000, s = 8, seed = 1)
dim(subs)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cate_demo_subsets", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cate_dnn")
### * cate_dnn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cate_dnn
### Title: CATE via DNN rank weights (orthogonalized local averaging)
### Aliases: cate_dnn

### ** Examples

set.seed(1)
n <- 50
x_mat <- matrix(rnorm(n * 2), n, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
psi <- rnorm(n)
out <- cate_dnn(
  X = x_mat, X0 = x0_mat, pseudo_outcome = psi,
  distance = "l2", method = "exact", k = 2L, s = 3L,
  se = FALSE, return = c("W", "ess")
)
out$tau



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cate_dnn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("confint.cate_result")
### * confint.cate_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: confint.cate_result
### Title: Confidence intervals for 'cate_result'
### Aliases: confint.cate_result

### ** Examples

set.seed(2)
n <- 25
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
confint(out)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("confint.cate_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dist_l1")
### * dist_l1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dist_l1
### Title: Fast L1 (Manhattan) distance to a single query point
### Aliases: dist_l1

### ** Examples

set.seed(1)
X <- matrix(rnorm(12), 4, 3)
x0 <- c(0, 1, -1)
dist_l1(X, x0)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dist_l1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dist_l1_scaled")
### * dist_l1_scaled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dist_l1_scaled
### Title: Scaled L1 distance (per-coordinate scaling)
### Aliases: dist_l1_scaled

### ** Examples

X <- matrix(rnorm(12), 4, 3)
dist_l1_scaled(X, c(0, 0, 0), inv_s = c(1, 2, 0.5))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dist_l1_scaled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dist_l2")
### * dist_l2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dist_l2
### Title: Fast L2 (Euclidean) distance to a single query point
### Aliases: dist_l2

### ** Examples

X <- matrix(rnorm(12), 4, 3)
dist_l2(X, c(0, 0, 0))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dist_l2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dist_l2_scaled")
### * dist_l2_scaled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dist_l2_scaled
### Title: Scaled L2 distance (per-coordinate scaling)
### Aliases: dist_l2_scaled

### ** Examples

X <- matrix(rnorm(12), 4, 3)
dist_l2_scaled(X, c(0, 0, 0), inv_s = c(1, 2, 0.5))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dist_l2_scaled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("distance_eval")
### * distance_eval

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: distance_eval
### Title: Evaluate a distance spec at a single target
### Aliases: distance_eval

### ** Examples

x_mat <- matrix(rnorm(12), 4, 3)
x0_vec <- c(0, 0, 0)
distance_eval(distance_l2(), x_mat, x0_vec)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("distance_eval", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("distance_eval_batch")
### * distance_eval_batch

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: distance_eval_batch
### Title: Evaluate a distance spec for many targets (batch)
### Aliases: distance_eval_batch

### ** Examples

x_mat <- matrix(rnorm(12), 4, 3)
x0_mat <- matrix(rnorm(6), 2, 3)
distance_eval_batch(distance_l2(), x_mat, x0_mat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("distance_eval_batch", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("distance_l1")
### * distance_l1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: distance_l1
### Title: L1 distance spec (single and batched evaluation)
### Aliases: distance_l1

### ** Examples

set.seed(3)
X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l1()
dspec$eval(X, X0[1, ])
dspec$eval_batch(X, X0)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("distance_l1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("distance_l1_scaled")
### * distance_l1_scaled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: distance_l1_scaled
### Title: Scaled L1 distance spec (one-time setup + eval/batch)
### Aliases: distance_l1_scaled

### ** Examples

X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l1_scaled(rule = "sd")
ctx <- dspec$setup(X)
dspec$eval_batch(X, X0, ctx)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("distance_l1_scaled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("distance_l2")
### * distance_l2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: distance_l2
### Title: L2 distance spec (single and batched evaluation)
### Aliases: distance_l2

### ** Examples

X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l2()
dspec$eval_batch(X, X0)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("distance_l2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("distance_l2_scaled")
### * distance_l2_scaled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: distance_l2_scaled
### Title: Scaled L2 distance spec (one-time setup + eval/batch)
### Aliases: distance_l2_scaled

### ** Examples

X <- matrix(rnorm(12), 4, 3)
X0 <- matrix(rnorm(6), 2, 3)
dspec <- distance_l2_scaled(rule = "sd")
ctx <- dspec$setup(X)
dspec$eval_batch(X, X0, ctx)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("distance_l2_scaled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("distance_names")
### * distance_names

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: distance_names
### Title: List registered distances
### Aliases: distance_names

### ** Examples

head(distance_names(), 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("distance_names", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dnn_weights_asymp_full_from_ranks")
### * dnn_weights_asymp_full_from_ranks

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dnn_weights_asymp_full_from_ranks", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dnn_weights_exact_full_from_ranks")
### * dnn_weights_exact_full_from_ranks

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dnn_weights_exact_full_from_ranks", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_distance")
### * get_distance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_distance
### Title: Retrieve a registered distance spec
### Aliases: get_distance

### ** Examples

d <- get_distance("l2")
is_distance(d)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_distance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_rank")
### * get_rank

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_rank
### Title: Retrieve a registered rank_weight
### Aliases: get_rank

### ** Examples

register_rank("tmp_rank_get", function(ranks) rep(1 / length(ranks), length(ranks)))
is_rank_weight(get_rank("tmp_rank_get"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_rank", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_score")
### * get_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_score
### Title: Retrieve a registered score
### Aliases: get_score

### ** Examples

sc0 <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "tmp_score_get")
register_score(sc0)
sc <- get_score("tmp_score_get")
is_score(sc)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is_compiled")
### * is_compiled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is_compiled
### Title: Does a score have a compiled implementation?
### Aliases: is_compiled

### ** Examples

is_compiled(new_score("r_only", order = 1L, eval = function(idxs, data, ctx, rng) data[idxs]))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is_compiled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is_distance")
### * is_distance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is_distance
### Title: Is this a distance spec?
### Aliases: is_distance

### ** Examples

is_distance(distance_l2())
is_distance(list(eval = function(X, x0, ctx) rep(0, nrow(X))))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is_distance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is_rank_weight")
### * is_rank_weight

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is_rank_weight
### Title: Is this a rank_weight spec?
### Aliases: is_rank_weight

### ** Examples

is_rank_weight(rank_weight_dnn_exact(k = 1, s = 1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is_rank_weight", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is_score")
### * is_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is_score
### Title: Is object a 'ustat_score'?
### Aliases: is_score

### ** Examples

sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
is_score(sc)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jackknife_fast_rank")
### * jackknife_fast_rank

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jackknife_fast_rank
### Title: Fast delete-one jackknife for rank-weighted functionals (auto
###   C++ when available, R fallback otherwise)
### Aliases: jackknife_fast_rank

### ** Examples

set.seed(1)
n <- 20
p <- 2
X <- matrix(rnorm(n * p), n, p)
ranks <- cbind(matrix(sample.int(n, n * 3, replace = TRUE), n, 3), sample.int(n, n, replace = TRUE))
a <- c(1, 1, 1)
phi <- cbind(rnorm(n))
out <- jackknife_fast_rank(ranks, a, phi, normalize = TRUE)
str(out$theta0)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jackknife_fast_rank", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jackknife_slow_rank")
### * jackknife_slow_rank

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jackknife_slow_rank
### Title: Slow delete-one jackknife for rank-weighted L-stat rules (fixed
###   ranks)
### Aliases: jackknife_slow_rank

### ** Examples

set.seed(1)
n <- 20
p <- 3
s <- 4
X <- matrix(rnorm(n * p), n, p)
# toy ranks: nearest-neighbor indices precomputed elsewhere (T x (s+1))
# here we fake something deterministic just for the example
ranks <- cbind(matrix(rep(1:s, each = n), n, s), s + 1L)
y <- rnorm(n)
a <- rep(1, s)
out <- jackknife_slow_rank(ranks, a, cbind(y))
str(out)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jackknife_slow_rank", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("list_distances")
### * list_distances

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: list_distances
### Title: Back-compat alias for 'distance_names()'
### Aliases: list_distances

### ** Examples

identical(distance_names(), list_distances())



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("list_distances", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("list_rank_weights")
### * list_rank_weights

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: list_rank_weights
### Title: Back-compat alias for 'rank_names()'
### Aliases: list_rank_weights

### ** Examples

identical(rank_names(), list_rank_weights())



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("list_rank_weights", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("list_ranks")
### * list_ranks

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: list_ranks
### Title: Back-compat alias for 'rank_names()'
### Aliases: list_ranks

### ** Examples

identical(rank_names(), list_ranks())



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("list_ranks", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("list_scores")
### * list_scores

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: list_scores
### Title: List registered scores
### Aliases: list_scores

### ** Examples

head(list_scores(), 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("list_scores", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("new_distance")
### * new_distance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: new_distance
### Title: Construct a distance spec
### Aliases: new_distance

### ** Examples

spec <- new_distance(
  name = "sq_l2",
  eval = function(X, x0, ctx = NULL) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
)
is_distance(spec)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("new_distance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("new_rank_weight")
### * new_rank_weight

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: new_rank_weight
### Title: Construct a rank-weight spec
### Aliases: new_rank_weight

### ** Examples

rw <- new_rank_weight(
  name = "inv_rank",
  eval = function(ranks, ctx = NULL) {
    w <- 1 / as.numeric(ranks)
    w / sum(w)
  }
)
is_rank_weight(rw)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("new_rank_weight", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("new_score")
### * new_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: new_score
### Title: Score specification (S3) for U-stat kernels
### Aliases: new_score

### ** Examples

sc <- new_score(
  name = "sum2",
  order = 2L,
  eval = function(idxs, data, ctx, rng) sum(data[idxs])
)
is_score(sc)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("new_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("np_demo_data")
### * np_demo_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: np_demo_data
### Title: Minimal nonparametric regression demo data (oracle m and Y only)
### Aliases: np_demo_data

### ** Examples

demo <- np_demo_data(n = 200, p = 2, design = "sinusoidal", seed = 1)
head(demo$data)
demo$m_x0




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("np_demo_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("np_demo_subsets")
### * np_demo_subsets

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: np_demo_subsets
### Title: Quick helper: sample a batch of subsets (singleton by default)
### Aliases: np_demo_subsets

### ** Examples

np_demo_subsets(n = 10, K = 4, s = 1, seed = 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("np_demo_subsets", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("np_dnn")
### * np_dnn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: np_dnn
### Title: Nonparametric regression via DNN rank weights
### Aliases: np_dnn

### ** Examples

set.seed(1)
n <- 60
x_mat <- matrix(rnorm(n * 2), n, 2)
y_vec <- x_mat[, 1] - 0.5 * x_mat[, 2] + rnorm(n, sd = 0.3)
out <- np_dnn(
  X = x_mat, x0 = c(0, 0), Y = y_vec,
  distance = "l2", method = "exact", k = 2L, s = 3L,
  se = FALSE, return = c("W", "ess")
)
out$mu



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("np_dnn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ord_to_ranks")
### * ord_to_ranks

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ord_to_ranks
### Title: Convert an index order to ranks (inverse permutation)
### Aliases: ord_to_ranks

### ** Examples

ord_to_ranks(c(3L, 1L, 2L))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ord_to_ranks", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("order_by_l1")
### * order_by_l1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: order_by_l1
### Title: Order rows by L1 distance to a point (stable ties)
### Aliases: order_by_l1

### ** Examples

set.seed(2)
X <- matrix(rnorm(10), 5, 2)
x0 <- c(0, 0)
order_by_l1(X, x0)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("order_by_l1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("order_by_l1_scaled")
### * order_by_l1_scaled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: order_by_l1_scaled
### Title: Order rows by scaled L1 (stable ties)
### Aliases: order_by_l1_scaled

### ** Examples

set.seed(6)
X <- matrix(rnorm(10), 5, 2)
x0 <- c(0, 0)
order_by_l1_scaled(X, x0, rule = "sd")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("order_by_l1_scaled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("order_by_l2")
### * order_by_l2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: order_by_l2
### Title: Order rows by L2 distance to a point (stable ties)
### Aliases: order_by_l2

### ** Examples

X <- matrix(rnorm(10), 5, 2)
order_by_l2(X, c(0, 0))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("order_by_l2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("order_by_l2_scaled")
### * order_by_l2_scaled

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: order_by_l2_scaled
### Title: Order rows by scaled L2 (stable ties)
### Aliases: order_by_l2_scaled

### ** Examples

X <- matrix(rnorm(10), 5, 2)
order_by_l2_scaled(X, c(0, 0), rule = "sd")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("order_by_l2_scaled", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plqr_demo_data")
### * plqr_demo_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plqr_demo_data
### Title: Minimal PLQR demo data (oracle nuisances; no internal fitting)
### Aliases: plqr_demo_data

### ** Examples

demo <- plqr_demo_data(
  n = 300, p = 2, p_d = 1, tau = 0.5,
  design = "sinusoidal", seed = 1
)
head(demo$data)
demo$truth$beta_true
demo$q_x0




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plqr_demo_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plqr_demo_subsets")
### * plqr_demo_subsets

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plqr_demo_subsets
### Title: Quick helper: sample a batch of subsets (for U-stat /
###   localization experiments)
### Aliases: plqr_demo_subsets

### ** Examples

plqr_demo_subsets(n = 10, K = 4, s = 3, seed = 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plqr_demo_subsets", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plqr_dnn")
### * plqr_dnn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plqr_dnn
### Title: PLQR local orthogonal moments via DNN rank weights
### Aliases: plqr_dnn

### ** Examples

set.seed(1)
demo <- plqr_demo_data(n = 120, p = 2, p_d = 2, tau = 0.5, seed = 1)
out <- plqr_dnn(
  data = demo$data,
  x_cols = c("X1", "X2"),
  D_cols = c("D1", "D2"),
  q_col = "q_oracle",
  mD_cols = c("mD1", "mD2"),
  y_col = "Y",
  beta = demo$truth$beta_true,
  tau = 0.5,
  x0 = c(0, 0),
  method = "exact",
  k = 2L, s = 3L,
  se = FALSE
)
out$moment



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plqr_dnn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plqr_moment")
### * plqr_moment

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plqr_moment
### Title: Orthogonal PLQR moment contributions for a given beta (oracle
###   nuisances)
### Aliases: plqr_moment

### ** Examples

demo <- plqr_demo_data(n = 100, p = 2, p_d = 2, tau = 0.5, seed = 1)
plqr_moment(demo$data, beta = demo$truth$beta_true, tau = 0.5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plqr_moment", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.cate_result")
### * predict.cate_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.cate_result
### Title: Predict from a 'cate_result' object at new targets
### Aliases: predict.cate_result

### ** Examples

set.seed(3)
n <- 40
x_mat <- matrix(rnorm(n * 2), n, 2)
fit <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
predict(fit, x0 = c(0.2, -0.1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.cate_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.cate_result")
### * print.cate_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.cate_result
### Title: Print a 'cate_result' object
### Aliases: print.cate_result

### ** Examples

set.seed(1)
n <- 20
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(
  X = x_mat,
  x0 = c(0, 0),
  pseudo_outcome = rnorm(n),
  se = FALSE
)
print(out)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.cate_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.score_diff_path")
### * print.score_diff_path

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.score_diff_path
### Title: Print a 'score_diff_path' object
### Aliases: print.score_diff_path

### ** Examples

sc <- new_score(
  name = "order1",
  order = 1L,
  eval_param = function(idxs, data, theta, ctx, rng) theta[1] * data[idxs]
)
print(score_diff_path(sc, data = c(1, 2, 3)))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.score_diff_path", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.ustat_score")
### * print.ustat_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.ustat_score
### Title: Print a 'ustat_score' object
### Aliases: print.ustat_score

### ** Examples

sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "sum2")
print(sc)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.ustat_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.weights_result")
### * print.weights_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.weights_result
### Title: Print a 'weights_result' object
### Aliases: print.weights_result

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
out <- weights_dnn_exact(x_mat, x0 = c(0, 0), k = 1L, s = 1L)
print(out)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.weights_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rank_names")
### * rank_names

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rank_names
### Title: List registered rank_weight names
### Aliases: rank_names

### ** Examples

head(rank_names(), 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rank_names", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rank_weight_dnn_asymp")
### * rank_weight_dnn_asymp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rank_weight_dnn_asymp
### Title: Rank-weight: DNN asymptotic (full U, WOR)
### Aliases: rank_weight_dnn_asymp

### ** Examples

rw <- rank_weight_dnn_asymp(k = 2, s = 3, normalize = TRUE)
rank_weight_eval(rw, c(2L, 1L, 3L, 4L))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rank_weight_dnn_asymp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rank_weight_dnn_exact")
### * rank_weight_dnn_exact

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rank_weight_dnn_exact
### Title: Rank-weight: DNN exact (full U, WOR)
### Aliases: rank_weight_dnn_exact

### ** Examples

rw <- rank_weight_dnn_exact(k = 2, s = 3)
rank_weight_eval(rw, c(2L, 1L, 3L, 4L))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rank_weight_dnn_exact", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rank_weight_eval")
### * rank_weight_eval

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rank_weight_eval
### Title: Evaluate rank-weight for a single target
### Aliases: rank_weight_eval

### ** Examples

rw <- as_rank_weight(function(r) rep(1 / length(r), length(r)))
rank_weight_eval(rw, c(2L, 1L, 3L))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rank_weight_eval", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rank_weight_eval_batch")
### * rank_weight_eval_batch

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rank_weight_eval_batch
### Title: Evaluate rank-weight for many targets (batch)
### Aliases: rank_weight_eval_batch

### ** Examples

rw <- as_rank_weight(function(r) rep(1 / length(r), length(r)))
rank_matrix <- cbind(c(1L, 2L, 3L), c(2L, 1L, 3L))
rank_weight_eval_batch(rw, rank_matrix)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rank_weight_eval_batch", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("register_distance")
### * register_distance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: register_distance
### Title: Register a distance spec under a name
### Aliases: register_distance

### ** Examples

nm <- "tmp_l1_demo"
register_distance(
  nm,
  function(X, x0, ctx = NULL) {
    rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  }
)
nm %in% distance_names()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("register_distance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("register_rank")
### * register_rank

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: register_rank
### Title: Register a rank_weight under a name
### Aliases: register_rank

### ** Examples

nm <- "tmp_rank_demo"
register_rank(nm, function(ranks) rep(1 / length(ranks), length(ranks)))
nm %in% rank_names()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("register_rank", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("register_score")
### * register_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: register_score
### Title: Register a score
### Aliases: register_score

### ** Examples

sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L, name = "tmp_score_demo")
register_score(sc)
"tmp_score_demo" %in% list_scores()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("register_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_cate_rho")
### * score_cate_rho

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_cate_rho
### Title: CATE score from a precomputed pseudo-outcome rho (pure score; no
###   weights)
### Aliases: score_cate_rho

### ** Examples

demo <- cate_demo_data(n = 200, p = 2, design = "hetero", seed = 1)
sc <- score_cate_rho(
  x0 = rep(0, 2), x_cols = paste0("X", 1:2),
  s = 5, k = 1, rho_col = "rho_oracle"
)
ctx <- score_setup(sc, demo$data)
subs <- cate_demo_subsets(n = nrow(demo$data), K = 3, s = 5, seed = 1)
score_eval_batch(sc, subs, demo$data, ctx) # three subset evaluations




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_cate_rho", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_diff_path")
### * score_diff_path

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_diff_path
### Title: Build a modular differentiable path for a score
### Aliases: score_diff_path

### ** Examples

sc <- new_score(
  name = "affine_order1",
  order = 1L,
  eval_param = function(idxs, data, theta, ctx, rng) {
    theta[1] * data[idxs] + theta[2]
  }
)
path <- score_diff_path(sc, data = c(1, 2, 4), method = "fd")
path$value(2L, theta = c(0.5, 1))
path$grad(2L, theta = c(0.5, 1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_diff_path", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_eval")
### * score_eval

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_eval
### Title: Evaluate kernel on one subset
### Aliases: score_eval

### ** Examples

sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
score_eval(sc, c(1L, 2L), c(3, 4, 5), NULL)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_eval", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_eval_batch")
### * score_eval_batch

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_eval_batch
### Title: Evaluate kernel on many subsets (R fast-path)
### Aliases: score_eval_batch

### ** Examples

sc <- as_score(function(idxs, data) sum(data[idxs]), s = 2L)
subs <- rbind(c(1L, 2L), c(2L, 3L))
score_eval_batch(sc, subs, c(3, 4, 5), NULL)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_eval_batch", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_names")
### * score_names

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_names
### Title: Back-compat alias for 'list_scores()'
### Aliases: score_names

### ** Examples

identical(score_names(), list_scores())



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_names", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_np_y")
### * score_np_y

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_np_y
### Title: Nonparametric regression *singleton* score from observed outcome
###   Y (order-1 pure score; no weights)
### Aliases: score_np_y

### ** Examples

set.seed(1)
n <- 50
p <- 2
dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
dat$Y <- with(dat, X1 - 2 * X2 + rnorm(n))

sc <- score_np_y(
  x0 = c(0, 0),
  x_cols = c("X1", "X2"),
  s = 1L, k = 1L,
  y_col = "Y"
)
ctx <- score_setup(sc, dat)
subs <- matrix(seq_len(n), ncol = 1L) # singleton subsets
score_eval_batch(sc, subs, dat, ctx) # returns dat$Y




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_np_y", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_order")
### * score_order

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_order
### Title: Resolve kernel order 's'
### Aliases: score_order

### ** Examples

sc <- new_score(
  name = "dyn",
  order = function(n, data, ctx) if (n > 5) 3L else 2L,
  eval = function(idxs, data, ctx, rng) sum(data[idxs])
)
score_order(sc, 4L, 1:4, NULL)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_order", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_plqr_component")
### * score_plqr_component

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_plqr_component
### Title: PLQR diagnostic score (component j); pure score (no weights, no
###   fitting)
### Aliases: score_plqr_component

### ** Examples

## Don't show: 
if (TRUE) withAutoprint({ # examplesIf
## End(Don't show)
demo <- plqr_demo_data(n = 200, p = 3, p_d = 2, tau = 0.5, seed = 1)
sc <- score_plqr_component(
  x0      = rep(0, 3),
  x_cols  = paste0("X", 1:3),
  D_cols  = paste0("D", 1:2),
  q_col   = "q_oracle",
  mD_cols = paste0("mD", 1:2),
  beta    = demo$truth$beta_true,
  j       = 1L,
  tau     = 0.5,
  s       = 5,
  k       = 1
)
ctx <- score_setup(sc, demo$data)
subs <- cate_demo_subsets(n = nrow(demo$data), K = 3, s = 5, seed = 1)
score_eval_batch(sc, subs, demo$data, ctx)
## Don't show: 
}) # examplesIf
## End(Don't show)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_plqr_component", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("score_setup")
### * score_setup

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: score_setup
### Title: Set up score context (precomputation)
### Aliases: score_setup

### ** Examples

sc <- new_score(
  name = "ctx_sum",
  order = 1L,
  setup = function(data, ...) list(offset = 1),
  eval = function(idxs, data, ctx, rng) data[idxs] + ctx$offset
)
ctx <- score_setup(sc, 1:3)
score_eval(sc, 2L, 1:3, ctx)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("score_setup", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_cate_data")
### * sim_cate_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_cate_data
### Title: Draw one sample from a named CATE simulation DGP
### Aliases: sim_cate_data

### ** Examples

dgp <- sim_cate_data("dgp2_nonlinear", n = 200, kx = 2, seed = 1)
head(dgp$data)
dgp$tau_x0




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_cate_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_cate_methods")
### * sim_cate_methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_cate_methods
### Title: List available CATE simulation methods
### Aliases: sim_cate_methods

### ** Examples

sim_cate_methods()




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_cate_methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_cate_one_rep")
### * sim_cate_one_rep

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_cate_one_rep
### Title: Run one replicate of CATE estimators
### Aliases: sim_cate_one_rep

### ** Examples

rep1 <- sim_cate_one_rep(
  n = 200, kx = 2, k = 5, s = 12,
  design = "dgp2_nonlinear",
  methods = c("dnn_oracle", "ols_interact"),
  seed = 1
)
rep1




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_cate_one_rep", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_cate_run")
### * sim_cate_run

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_cate_run
### Title: Run the CATE simulation grid
### Aliases: sim_cate_run

### ** Examples

## Not run: 
##D # Single method (default):
##D raw <- sim_cate_run(designs = "dgp2_nonlinear", n_vals = 500L, n_reps = 10L)
##D 
##D # Multiple methods on the same draws:
##D raw <- sim_cate_run(
##D   designs = "dgp2_nonlinear", n_vals = 500L, n_reps = 10L,
##D   methods = c("dnn_oracle", "dnn_crossfit", "ols_interact")
##D )
##D summarise_sim_results(raw)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_cate_run", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summarise_sim_results")
### * summarise_sim_results

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summarise_sim_results
### Title: Summarise raw per-replicate CATE simulation results
### Aliases: summarise_sim_results

### ** Examples

## Not run: 
##D raw <- sim_cate_run(designs = "dgp1_linear", n_vals = 500L, n_reps = 10L)
##D summarise_sim_results(raw)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summarise_sim_results", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.cate_result")
### * summary.cate_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.cate_result
### Title: Summarize a 'cate_result' object
### Aliases: summary.cate_result

### ** Examples

set.seed(1)
n <- 30
x_mat <- matrix(rnorm(n * 2), n, 2)
out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n))
summary(out)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.cate_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("weights_dnn")
### * weights_dnn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: weights_dnn
### Title: DNN weights (choose method), distance-based ranking
### Aliases: weights_dnn

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0 <- c(0, 0)
out <- weights_dnn(x_mat, x0 = x0, k = 2L, s = 3L, method = "exact")
dim(out$W)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("weights_dnn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("weights_dnn_asymp")
### * weights_dnn_asymp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: weights_dnn_asymp
### Title: DNN weights (asymptotic full-U), distance-based ranking
### Aliases: weights_dnn_asymp

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(4), 2, 2)
out <- weights_dnn_asymp(x_mat, x0_mat, k = 2L, s = 3L, normalize = TRUE)
dim(out$W)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("weights_dnn_asymp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("weights_dnn_exact")
### * weights_dnn_exact

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: weights_dnn_exact
### Title: DNN weights (exact full-U), distance-based ranking
### Aliases: weights_dnn_exact

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
out <- weights_dnn_exact(x_mat, x0_mat, distance = "l2", k = 2L, s = 3L, return = "W")
dim(out$W)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("weights_dnn_exact", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("weights_matrix")
### * weights_matrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: weights_matrix
### Title: Compute weight matrix W via an arbitrary ranking + weight
###   mapping
### Aliases: weights_matrix

### ** Examples

set.seed(1)
x_mat <- matrix(rnorm(20), 10, 2)
x0_mat <- matrix(rnorm(6), 3, 2)
rw <- as_rank_weight(function(r) {
  w <- 1 / as.numeric(r)
  w / sum(w)
}, name = "inv_rank")
out <- weights_matrix(x_mat, x0_mat, ranking = "l2", weight = rw, return = c("W", "ranks"))
dim(out$W)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("weights_matrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
