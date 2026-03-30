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




