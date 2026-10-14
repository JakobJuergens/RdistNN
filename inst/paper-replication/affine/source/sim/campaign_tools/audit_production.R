#!/usr/bin/env Rscript
# Audit completed production files only. No data or estimator is regenerated.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE))
code <- normalizePath(file.path(dirname(script),".."))
source(file.path(code,"../bootstrap.R"))
for(f in c("affine_config.R","affine_models.R","affine_series.R","affine_local.R","affine_kernel.R",
           "affine_diagnostics.R","affine_summary.R","affine_runner.R")) source(file.path(code,f))
args <- commandArgs(TRUE)
if(length(args)!=2L) stop("Usage: audit_production.R CAMPAIGN_ROOT NEW_EXPORT_DIR")
campaign <- normalizePath(args[1]); out <- args[2]
if(file.exists(out)) stop("Choose a new audit export directory")
proposal <- jsonlite::read_json(file.path(campaign,"configs/proposal.json"),simplifyVector=TRUE)
stopifnot(all(vapply(proposal$profiles,function(p) p$replications==10000,logical(1))))
# Refuse a partial campaign before reading or exporting large payloads.
for(name in names(proposal$profiles)) {
  m <- jsonlite::read_json(file.path(campaign,"runs",name,"manifest.json"),simplifyVector=TRUE)
  attempts <- list.dirs(file.path(campaign,"runs/measurements"),recursive=FALSE,full.names=TRUE)
  attempts <- attempts[startsWith(basename(attempts),paste0(name,"_"))]
  processes <- lapply(attempts,function(a) jsonlite::read_json(file.path(a,"process.json"),simplifyVector=TRUE))
  stopifnot(length(processes)>0)
  p <- processes[[which.max(vapply(processes,`[[`,numeric(1),"started_unix"))]]
  stopifnot(m$status=="completed",m$completed_replications==10000,p$status=="completed",p$returncode==0)
}
# Match the recorded override string as well as the actual source files. The
# runbook uses an absolute path; an equivalent relative path changes this
# metadata field and must not cause a spurious runtime-freeze mismatch.
first <- readRDS(file.path(campaign,"runs",names(proposal$profiles)[1],"freeze.rds"))
override <- first$details$dependency$source_override
if(nzchar(override)) Sys.setenv(SIM_RDISTNN_SOURCE_DIR=override) else Sys.unsetenv("SIM_RDISTNN_SOURCE_DIR")
affine_load_dependency();dir.create(out,recursive=TRUE)
receipts <- tables <- list(); shared_seeds <- NULL; overlap <- list()
eq <- function(a,b,tol=1e-10) stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))
for(name in names(proposal$profiles)) {
  path <- file.path(campaign,"runs",name);cfg <- readRDS(file.path(path,"config.rds"))
  spec <- proposal$profiles[[name]]; freeze <- readRDS(file.path(path,"freeze.rds"))
  manifest <- jsonlite::read_json(file.path(path,"manifest.json"),simplifyVector=TRUE)
  schema <- jsonlite::read_json(file.path(path,"schema.json"),simplifyVector=TRUE)
  seeds <- affine_seed_map(cfg$seed_capacity,cfg$master_seed,cfg$stage)
  stopifnot(cfg$stage=="production",all(cfg$cells$replications==10000),cfg$replications==10000,
    identical(freeze$hash,affine_freeze(cfg,code,seeds)$hash),isTRUE(schema$validation$raw_fields_valid),
    identical(digest::digest(file=file.path(path,"config.rds"),algo="sha256"),proposal$config_hashes[[name]]))
  if(is.null(shared_seeds)) shared_seeds <- seeds else stopifnot(identical(seeds,shared_seeds))
  pilot <- affine_seed_map(cfg$seed_capacity,cfg$master_seed,"pilot")
  stopifnot(!length(intersect(unlist(seeds[-1]),unlist(pilot[-1]))))
  raw <- readRDS(file.path(path,"raw.rds"));summary <- readRDS(file.path(path,"summary.rds"))
  stopifnot(nrow(raw)==spec$expected_main_rows+spec$expected_control_rows,
    nrow(raw)==manifest$actual_rows,!anyDuplicated(affine_row_keys(raw)),
    identical(digest::digest(raw,algo="sha256"),manifest$numerical_hash),
    identical(digest::digest(summary,algo="sha256"),manifest$summary_hash))
  expected <- affine_expected_keys(cfg,1:10000)
  stopifnot(setequal(affine_row_keys(raw),affine_row_keys(expected)))
  # Stream sealed chunks; verify hashes, IDs and exact raw correspondence without
  # retaining a second copy of all diagnostic families. Native field/ESS checks
  # already passed in the completed runner's full exported-schema validation.
  chunk_ids <- integer(); chunk_count <- 0L
  raw_keys <- affine_row_keys(raw)
  for(file in list.files(file.path(path,"chunks"),pattern="\\.rds$",full.names=TRUE)) {
    ch <- readRDS(file);chunk_count <- chunk_count+1L
    stopifnot(identical(ch$freeze_hash,freeze$hash),identical(ch$payload_hash,affine_payload_hash(ch)),
      identical(ch$raw,do.call(rbind,lapply(ch$families,`[[`,"raw"))))
    ix <- match(affine_row_keys(ch$raw),raw_keys);stopifnot(!anyNA(ix))
    a <- ch$raw;b <- raw[ix,];rownames(a) <- rownames(b) <- NULL
    stopifnot(identical(a,b));chunk_ids <- c(chunk_ids,ch$replication_ids)
    rm(ch,a,b)
  }
  stopifnot(!anyDuplicated(chunk_ids),identical(sort(chunk_ids),1:10000),chunk_count==1250)
  for(i in seq_len(nrow(summary$performance))) {
    p <- summary$performance[i,]
    z <- raw[raw$cell_id==p$cell_id & raw$application==p$application &
      raw$training_method==p$training_method & raw$localization_method==p$localization_method,]
    good <- z$status=="ok" & is.finite(z$estimate)
    route <- p$variance_method;se <- z[[paste0("SE_",route)]]
    lo <- z[[paste0("lower_",route)]];hi <- z[[paste0("upper_",route)]]
    valid <- good & z[[paste0("variance_status_",route)]]!="unavailable" &
      is.finite(se) & se>=0 & is.finite(lo) & is.finite(hi)
    stopifnot(p$successful_estimates==sum(good),p$intervals==sum(valid))
    eq(p$bias,if(any(good)) mean(z$estimate[good]-z$point_target[good]) else NA_real_)
    eq(p$sd,if(sum(good)>1) sd(z$estimate[good]) else NA_real_)
    eq(p$coverage,if(any(valid)) mean(lo[valid]<=z$point_target[valid] & hi[valid]>=z$point_target[valid]) else NA_real_)
    eq(p$mean_length,if(any(valid)) mean(hi[valid]-lo[valid]) else NA_real_)
  }
  dg <- summary$diagnostics
  maximum <- function(v) if(any(is.finite(v))) max(abs(v[is.finite(v)])) else NA_real_
  identity <- maximum(dg$decomposition_identity_error)
  covariance <- maximum(dg$covariance_transfer-dg$covariance_component-dg$covariance_target-dg$covariance_cross)
  stopifnot(is.na(identity)||identity<1e-10,is.na(covariance)||covariance<1e-10)
  for(k in names(summary)[vapply(summary,is.data.frame,logical(1))]) {
    z <- summary[[k]]
    if(k=="diagnostics") next # full diagnostic tables remain in local run files
    if(nrow(z) && "cell_id"%in%names(z)) {
      cells <- cfg$cells[match(z$cell_id,cfg$cells$cell_id),]
      z <- cbind(profile=name,z,cells[c("n","design_id","order_mult","M","K")])
    }
    tables[[k]][[name]] <- z
  }
  routes <- c("jackknife","direct_moment","direct_squared")
  receipts[[name]] <- list(rows=nrow(raw),main_repetitions=10000,cells=nrow(cfg$cells),
    validated_chunks=chunk_count,point_failures=sum(raw$status!="ok"),
    diagnostic_failures=sum(raw$diagnostic_status!="ok"),
    unavailable_variances=as.list(setNames(vapply(routes,function(v)sum(raw[[paste0("variance_status_",v)]]=="unavailable"),integer(1)),routes)),
    largest_equation_residual=maximum(raw$score_mean_residual),largest_gap_identity_error=identity,
    largest_covariance_identity_error=covariance,freeze_hash=freeze$hash,
    raw_hash=manifest$numerical_hash,summary_hash=manifest$summary_hash,worker_cpu_seconds=manifest$worker_cpu_seconds)
  if(name%in%c("dnn","kernel_comparison")) {
    z <- raw[raw$training_method%in%c("full","pooled_five_fold","oracle") & raw$localization_method=="dnn",]
    fields <- c("cell_id","replication_id","application","training_method","localization_method","data_seed","fold_seed",
      "estimate","Gamma_hat","N",grep("^(variance|SE|lower|upper|score_cov_scaled)_",names(z),value=TRUE))
    overlap[[name]] <- z[fields]
  }
  cat(name,": audited",nrow(raw),"rows and",chunk_count,"chunks\n")
  rm(raw,summary,dg,expected);gc()
}
a <- overlap$dnn;b <- overlap$kernel_comparison
pos <- match(affine_row_keys(b),affine_row_keys(a));keep <- !is.na(pos)
stopifnot(sum(keep)==480000)
for(k in setdiff(names(a),c("cell_id","replication_id","application","training_method","localization_method")))
  stopifnot(identical(a[[k]][pos[keep]],b[[k]][keep]))
for(k in names(tables)) {
  z <- tables[[k]][vapply(tables[[k]],nrow,integer(1))>0]
  if(length(z)) write.csv(do.call(rbind,z),file.path(out,paste0(k,".csv")),row.names=FALSE)
}
receipt <- list(status="completed_and_audited_production",total_rows=sum(vapply(receipts,`[[`,numeric(1),"rows")),
  profiles=receipts,common_independent_draw_families=10000,duplicate_DNN_rows=480000,
  duplicate_DNN_estimates_and_variances_identical=TRUE,production_pilot_banks_disjoint=TRUE,
  all_planned_keys_and_chunk_raw_correspondence_verified=TRUE,
  audit_script_sha256=digest::digest(file=script,algo="sha256"),
  manuscript_promotion="not_performed")
stopifnot(receipt$total_rows==proposal$expected_total_rows)
jsonlite::write_json(receipt,file.path(out,"audit.json"),pretty=TRUE,auto_unbox=TRUE,digits=17,na="null")
cat("Full 10000-repetition campaign audited; no estimation or manuscript promotion performed.\n")
