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



