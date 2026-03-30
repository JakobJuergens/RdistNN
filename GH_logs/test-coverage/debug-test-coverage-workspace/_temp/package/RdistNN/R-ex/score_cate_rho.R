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




