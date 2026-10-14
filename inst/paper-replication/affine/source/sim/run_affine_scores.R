#!/usr/bin/env Rscript
# Standalone entry point; never writes into a historical campaign directory.
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
if(length(script_arg)!=1L) stop("Run this entry point with Rscript")
affine_source_dir <- dirname(normalizePath(sub("^--file=","",script_arg)))
source(file.path(dirname(affine_source_dir),"bootstrap.R"))
for(file in c("affine_config.R","affine_models.R","affine_series.R","affine_local.R","affine_kernel.R",
             "affine_diagnostics.R","affine_summary.R","affine_runner.R")) source(file.path(affine_source_dir,file))
if(!nzchar(Sys.getenv("TZ",unset=""))) Sys.setenv(TZ="UTC")
config_file <- Sys.getenv("SIM_AFFINE_CONFIG_RDS",unset="")
cfg <- if(nzchar(config_file)) readRDS(config_file) else affine_config(Sys.getenv("SIM_AFFINE_STAGE",unset="preflight"),
    list(study=Sys.getenv("SIM_AFFINE_STUDY",unset="dnn")))
number_env <- function(name,default,integer=FALSE) {
  value <- Sys.getenv(name,unset="")
  if(!nzchar(value)) return(default)
  out <- suppressWarnings(as.numeric(value))
  if(length(out)!=1L || !is.finite(out) || out<=0 || (integer && out!=floor(out))) stop(paste("Invalid",name))
  if(integer) as.integer(out) else out
}
cfg$workers <- number_env("SIM_AFFINE_WORKERS",cfg$workers,TRUE)
cfg$chunk_size <- number_env("SIM_AFFINE_CHUNK_SIZE",cfg$chunk_size,TRUE)
cfg$cpu_budget <- number_env("SIM_AFFINE_CPU_BUDGET",cfg$cpu_budget)
if(nzchar(Sys.getenv("SIM_AFFINE_REPLICATIONS",unset=""))) {
  replications <- number_env("SIM_AFFINE_REPLICATIONS",cfg$replications,TRUE)
  if(replications>cfg$seed_capacity) stop("Replications exceed the frozen seed capacity")
  cfg$replications <- replications; cfg$cells$replications <- replications
  cfg$cells$control_replications <- pmin(cfg$control_reps,replications)
}
run_name <- Sys.getenv("SIM_AFFINE_RUN_NAME",unset="")
if(nzchar(run_name)) {
  if(!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$",run_name)) stop("Run name must be a simple directory identifier")
  cfg$campaign_id <- run_name
} else run_name <- cfg$campaign_id
output_dir <- Sys.getenv("SIM_AFFINE_OUTPUT_DIR",unset="")
if(!nzchar(output_dir)) output_dir <- file.path(dirname(sim_find_project_root(affine_source_dir)),
       "02_code","simulations_cond_z_est",run_name)
resume_text <- tolower(Sys.getenv("SIM_AFFINE_RESUME",unset="false"))
if(!resume_text%in%c("true","false","1","0","yes","no")) stop("Invalid SIM_AFFINE_RESUME")
max_chunks <- number_env("SIM_AFFINE_MAX_CHUNKS",Inf,TRUE)
result <- affine_run(cfg,output_dir,resume=resume_text%in%c("true","1","yes"),
       max_chunks=max_chunks,source_dir=affine_source_dir)
cat(sprintf("Affine %s: %s; %d rows; %d numerical failures; %.3f worker CPU seconds\n",
    cfg$stage,result$status,nrow(result$raw),result$manifest$numerical_failures,result$manifest$worker_cpu_seconds))
cat("Receipt:",normalizePath(output_dir),"\n")
