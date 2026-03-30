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




