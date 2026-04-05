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




