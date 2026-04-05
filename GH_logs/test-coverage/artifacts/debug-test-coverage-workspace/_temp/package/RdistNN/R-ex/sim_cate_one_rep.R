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




