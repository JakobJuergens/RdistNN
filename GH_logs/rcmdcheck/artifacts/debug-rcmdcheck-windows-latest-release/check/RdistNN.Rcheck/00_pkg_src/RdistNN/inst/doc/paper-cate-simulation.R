## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)
library(RdistNN)

## ----replication-path, eval = FALSE-------------------------------------------
# system.file("paper-replication/cate/run.R", package = "RdistNN")

## ----dgp-example--------------------------------------------------------------
dgp <- sim_cate_data("dgp2_nonlinear", n = 300, kx = 2, seed = 42)
head(dgp$data[, c("X1", "X2", "D", "Y", "tau_oracle", "e_oracle")])
dgp$tau_x0   # true CATE at x0 = (0, 0)

## ----one-rep------------------------------------------------------------------
r <- sim_cate_one_rep(
  n      = 300,
  kx     = 2,
  k      = 5,
  s      = floor(300 ^ (3 / 4)),   # undersmoothing: s ~ n^((kx+1)/(kx+2))
  design = "dgp2_nonlinear",
  seed   = 1
)
r

## ----small-run----------------------------------------------------------------
raw <- sim_cate_run(
  designs     = c("dgp1_linear", "dgp2_nonlinear"),
  n_vals      = 300L,
  kx_vals     = 2L,
  n_reps      = 5L,
  master_seed = 42L,
  verbose     = FALSE
)
head(raw)

## ----summary------------------------------------------------------------------
tab <- summarise_sim_results(raw)
print(tab, digits = 3)

## ----parallel-pc, eval = FALSE------------------------------------------------
# library(future)
# library(future.apply)
# 
# plan(multisession)            # use all local cores
# raw <- sim_cate_run(
#   designs = "dgp2_nonlinear",
#   n_vals  = 500L,
#   n_reps  = 50L
# )
# plan(sequential)              # reset when done

## ----hpc-path, eval = FALSE---------------------------------------------------
# # Locate the HPC scripts after installing the package:
# system.file("paper-replication/cate/hpc", package = "RdistNN")

## ----full-run, eval = FALSE---------------------------------------------------
# source(system.file("paper-replication/cate/run.R", package = "RdistNN"))

## ----benchmarks-path, eval = FALSE--------------------------------------------
# system.file("paper-replication/cate/benchmarks.R", package = "RdistNN")

