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



