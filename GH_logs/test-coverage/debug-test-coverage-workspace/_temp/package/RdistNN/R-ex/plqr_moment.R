### Name: plqr_moment
### Title: Orthogonal PLQR moment contributions for a given beta (oracle
###   nuisances)
### Aliases: plqr_moment

### ** Examples

demo <- plqr_demo_data(n = 100, p = 2, p_d = 2, tau = 0.5, seed = 1)
plqr_moment(demo$data, beta = demo$truth$beta_true, tau = 0.5)



