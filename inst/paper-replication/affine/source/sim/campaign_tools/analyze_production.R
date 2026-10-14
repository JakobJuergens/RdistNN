#!/usr/bin/env Rscript
# Analysis of saved production only; never generates data or fits an estimator.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
code <- normalizePath(file.path(dirname(script), ".."))
source(file.path(code, "../bootstrap.R"))
for(f in c("affine_config.R","affine_models.R","affine_series.R","affine_local.R",
           "affine_kernel.R","affine_diagnostics.R","affine_summary.R","affine_runner.R"))
  source(file.path(code,f))
args <- commandArgs(TRUE)
if(length(args)!=3L) stop("Usage: analyze_production.R CAMPAIGN_ROOT AUDITED_EXPORT ANALYSIS_DATA_DIR")
campaign <- normalizePath(args[1]); accepted <- normalizePath(args[2]); out <- args[3]
audit <- jsonlite::read_json(file.path(accepted,"audit.json"),simplifyVector=TRUE)
stopifnot(audit$status=="completed_and_audited_production",audit$total_rows==2516400)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
methods <- c("full","pooled_five_fold","oracle")
routes <- c("jackknife","direct_moment","direct_squared")
keys <- c("cell_id","application","training_method","localization_method")
probs <- sort(unique(c(.001,.0025,.005,.01,.025,seq(.05,.95,.01),.975,.99,.995,.9975,.999)))
moments <- function(v) {
  v <- v[is.finite(v)]; n <- length(v)
  if(!n) return(data.frame(observed=0,mean=NA_real_,sd=NA_real_,mcse=NA_real_,median=NA_real_,q05=NA_real_,q95=NA_real_,maximum=NA_real_))
  data.frame(observed=n,mean=mean(v),sd=if(n>1) sd(v) else NA_real_,
    mcse=if(n>1) sd(v)/sqrt(n) else NA_real_,median=affine_quantile(v,.5),
    q05=affine_quantile(v,.05),q95=affine_quantile(v,.95),maximum=max(v))
}
shape <- function(v) {
  v <- v[is.finite(v)]; a <- v-mean(v); m2 <- mean(a^2)
  c(mean=mean(v),sd=sd(v),skewness=if(m2>0) mean(a^3)/m2^1.5 else NA_real_,
    excess_kurtosis=if(m2>0) mean(a^4)/m2^2-3 else NA_real_)
}
metadata <- function(z) z[1,c("profile",keys,"design_id","n","M","K","s","lambda","order_mult","canonical"),drop=FALSE]
metrics <- qq <- dgstats <- gram <- banks <- receipts <- list()
dnn_cells <- readRDS(file.path(campaign,"runs/dnn/config.rds"))$cells$cell_id
for(name in names(audit$profiles)) {
  path <- file.path(campaign,"runs",name); cfg <- readRDS(file.path(path,"config.rds"))
  manifest <- jsonlite::read_json(file.path(path,"manifest.json"),simplifyVector=TRUE)
  raw <- readRDS(file.path(path,"raw.rds")); summary <- readRDS(file.path(path,"summary.rds"))
  stopifnot(manifest$status=="completed",manifest$completed_replications==10000,
    identical(digest::digest(raw,algo="sha256"),manifest$numerical_hash),
    identical(digest::digest(summary,algo="sha256"),manifest$summary_hash),
    identical(manifest$numerical_hash,audit$profiles[[name]]$raw_hash),
    identical(manifest$summary_hash,audit$profiles[[name]]$summary_hash))
  dg <- summary$diagnostics
  stopifnot(identical(affine_row_keys(raw),affine_row_keys(dg)))
  raw$profile <- name; raw$order_mult <- cfg$cells$order_mult[match(raw$cell_id,cfg$cells$cell_id)]
  raw$canonical <- !(name=="kernel_comparison" & raw$localization_method=="dnn" & raw$cell_id%in%dnn_cells)
  dg$covariance_oracle_cross <- dg$fitted_score_scaled-dg$covariance_point_oracle_scaled-dg$covariance_transfer
  for(f in c("transfer","component","target","cross","oracle_cross")) {
    den <- dg$covariance_point_oracle_scaled
    dg[[paste0("covariance_",f,"_ratio_point_oracle")]] <- ifelse(is.finite(den)&den>0,dg[[paste0("covariance_",f)]]/den,NA_real_)
  }
  dg$point_cate_prediction_MSE_sum <- dg$fresh_pi+dg$fresh_mu0+dg$fresh_mu1
  dg$point_slope_prediction_MSE_sum <- dg$fresh_pi+dg$fresh_m
  dg$point_cate_product_bound <- sqrt(dg$fresh_pi)*(sqrt(dg$fresh_mu0)+sqrt(dg$fresh_mu1))
  dg$point_slope_product_bound <- abs(1/8)*dg$fresh_pi+sqrt(dg$fresh_pi*dg$fresh_m)
  groups <- split(seq_len(nrow(raw)),interaction(raw[keys],drop=TRUE,lex.order=TRUE))
  for(ix in groups) {
    z <- raw[ix,,drop=FALSE]; d <- dg[ix,,drop=FALSE]; meta <- metadata(z)
    good <- z$status=="ok" & is.finite(z$estimate)
    mean_est <- mean(z$estimate[good]); mc_sd <- sd(z$estimate[good])
    for(route in routes) {
      se <- z[[paste0("SE_",route)]]; lo <- z[[paste0("lower_",route)]]; hi <- z[[paste0("upper_",route)]]
      valid <- good & z[[paste0("variance_status_",route)]]!="unavailable" & is.finite(se)&se>0&is.finite(lo)&is.finite(hi)
      valid[is.na(valid)] <- FALSE
      point_t <- (z$estimate[valid]-z$point_target[valid])/se[valid]
      centered_t <- (z$estimate[valid]-mean_est)/se[valid]
      sd_z <- (z$estimate[valid]-mean_est)/mc_sd
      metrics[[length(metrics)+1L]] <- cbind(meta,data.frame(variance_method=route,completed=nrow(z),
        successful=sum(good),studentized=sum(valid),excluded=nrow(z)-sum(valid),
        mean_N=mean(z$N[good]),min_N=min(z$N[good]),max_N=max(z$N[good]),
        mean_Gamma=mean(z$Gamma_hat[good]),Gamma_mcse=sd(z$Gamma_hat[good])/sqrt(sum(good)),
        mean_Gamma_minus_point=mean(z$Gamma_hat[good])-if(z$application[1]=="cate") 1 else 31/128,
        mean_radius=mean(z$weighted_radius[good]),radius_q95=affine_quantile(z$weighted_radius[good],.95),
        mean_largest_weight=mean(z$largest_weight[good]),bias=mean(z$estimate[good]-z$point_target[good]),
        bias_mcse=mc_sd/sqrt(sum(good)),sd=mc_sd,
        coverage_point=mean(abs(point_t)<=qnorm(.975)),
        coverage_empirical_center=mean(abs(centered_t)<=qnorm(.975)),
        coverage_realized_weighted_truth=mean(lo[valid]<=z$weighted_truth[valid]&hi[valid]>=z$weighted_truth[valid]),
        as.list(setNames(shape(point_t),paste0("point_t_",names(shape(point_t))))),
        as.list(setNames(shape(centered_t),paste0("centered_t_",names(shape(centered_t))))),
        as.list(setNames(shape(sd_z),paste0("mc_z_",names(shape(sd_z)))))))
      if(z$canonical[1] && z$training_method[1]%in%methods) for(definition in c("point_t","centered_t","mc_z")) {
        value <- switch(definition,point_t=point_t,centered_t=centered_t,mc_z=sd_z)
        qq[[length(qq)+1L]] <- cbind(meta[rep(1,length(probs)),],data.frame(variance_method=route,
          definition=definition,observed=length(value),probability=probs,normal_quantile=qnorm(probs),
          empirical_quantile=affine_quantile(value,probs)))
      }
    }
    if(z$canonical[1]) for(field in names(d)[vapply(d,is.numeric,logical(1))]) {
      dgstats[[length(dgstats)+1L]] <- cbind(meta,data.frame(field=field),moments(d[[field]][good]))
    }
    if(z$canonical[1] && z$training_method[1]=="full" && z$localization_method[1]=="dnn" && any(!is.na(z$gram_event))) {
      for(event in c(TRUE,FALSE)) {
        v <- good & !is.na(z$gram_event)&z$gram_event==event
        if(!any(v)) next
        cover <- z$lower_jackknife[v]<=z$point_target[v]&z$upper_jackknife[v]>=z$point_target[v]
        ci <- affine_binomial_interval(sum(cover),length(cover))
        gram[[length(gram)+1L]] <- cbind(meta,data.frame(gram_event=event,observed=sum(v),
          bias=mean(z$estimate[v]-z$point_target[v]),bias_mcse=sd(z$estimate[v]-z$point_target[v])/sqrt(sum(v)),
          coverage=mean(cover),coverage_lower=ci[1],coverage_upper=ci[2],
          mean_Gamma=mean(z$Gamma_hat[v]),min_raw_gram=min(z$min_raw_gram[v]),
          mean_clip=mean(z$max_clip_fraction[v])))
      }
    }
  }
  fields <- unique(c("profile","canonical",keys,"replication_id","data_seed","fold_seed","integration_seed",
    "design_id","n","M","K","s","lambda","order_mult","point_target","estimate","status","N","N_reference",
    grep("^(SE|lower|upper|variance_status)_",names(raw),value=TRUE)))
  banks[[name]] <- raw[raw$canonical & raw$training_method%in%methods,fields]
  receipts[[name]] <- list(raw_hash=manifest$numerical_hash,summary_hash=manifest$summary_hash,
    raw_file_sha256=digest::digest(file=file.path(path,"raw.rds"),algo="sha256"),
    summary_file_sha256=digest::digest(file=file.path(path,"summary.rds"),algo="sha256"),rows=nrow(raw))
  cat(name,": extracted saved raw/diagnostic summaries\n")
  rm(raw,summary,dg,groups);gc()
}
bank <- do.call(rbind,banks);rm(banks);gc()
pairs <- list()
compare <- function(a,b,kind,left_route="jackknife",right_route=left_route) {
  if(!nrow(a)||!nrow(b)) return(invisible(NULL))
  stopifnot(!anyDuplicated(a$replication_id),!anyDuplicated(b$replication_id))
  b <- b[match(a$replication_id,b$replication_id),];stopifnot(!anyNA(b$replication_id))
  for(k in c("data_seed","fold_seed","integration_seed","point_target")) stopifnot(identical(a[[k]],b[[k]]))
  good <- a$status=="ok" & b$status=="ok" & is.finite(a$estimate)&is.finite(b$estimate)
  available <- function(z,r) z[[paste0("variance_status_",r)]]!="unavailable" &
    is.finite(z[[paste0("lower_",r)]])&is.finite(z[[paste0("upper_",r)]])
  valid <- good & available(a,left_route)&available(b,right_route);valid[is.na(valid)] <- FALSE
  coverage <- function(z,r) as.integer(z[[paste0("lower_",r)]][valid]<=z$point_target[valid]&z[[paste0("upper_",r)]][valid]>=z$point_target[valid])
  length <- function(z,r) z[[paste0("upper_",r)]][valid]-z[[paste0("lower_",r)]][valid]
  gap <- a$estimate[good]-b$estimate[good]
  loss <- (a$estimate[good]-a$point_target[good])^2-(b$estimate[good]-b$point_target[good])^2
  mcse <- function(v) if(base::length(v)>1) sd(v)/sqrt(base::length(v)) else NA_real_
  pairs[[base::length(pairs)+1L]] <<- cbind(a[1,c("design_id","n","M","K","s","order_mult","application","training_method","localization_method")],
    data.frame(comparison=kind,left_cell=a$cell_id[1],right_cell=b$cell_id[1],right_design=b$design_id[1],right_n=b$n[1],
      right_M=b$M[1],right_training=b$training_method[1],right_localizer=b$localization_method[1],
      left_route=left_route,right_route=right_route,point_pairs=sum(good),excluded_points=nrow(a)-sum(good),
      mean_gap=mean(gap),mean_gap_mcse=mcse(gap),RMS_gap=sqrt(mean(gap^2)),
      MSE_difference=mean(loss),MSE_difference_mcse=mcse(loss),paired=sum(valid),excluded=nrow(a)-sum(valid),
      left_coverage=mean(coverage(a,left_route)),right_coverage=mean(coverage(b,right_route))),
    affine_coverage_comparison(coverage(a,left_route)-coverage(b,right_route),length(a,left_route)-length(b,right_route)))
}
groups <- split(bank,interaction(bank[c("design_id","n","order_mult","application","training_method","localization_method")],drop=TRUE))
for(a in groups) {
  d <- a$design_id[1];n <- a$n[1];o <- a$order_mult[1];app <- a$application[1];t <- a$training_method[1];loc <- a$localization_method[1]
  pick <- function(design=d,size=n,order=o,training=t,localizer=loc)
    bank[bank$design_id==design & bank$n==size & bank$order_mult==order & bank$application==app &
      bank$training_method==training & bank$localization_method==localizer,]
  if(t=="full") for(t2 in c("pooled_five_fold","oracle")) compare(a,pick(training=t2),paste0("full_minus_",t2))
  for(r in c("direct_moment","direct_squared")) compare(a,a,paste0("jackknife_minus_",r),"jackknife",r)
  if(loc=="dnn") {
    b <- pick(localizer="kernel_gaussian")
    compare(a,b,"dnn_jk_minus_kernel_jk")
    compare(a,b,"dnn_jk_minus_kernel_sandwich","jackknife","direct_squared")
  }
  if(d%in%c("C0H","C2S")) {
    b <- pick(design=if(d=="C0H") "C0" else "C2")
    compare(a,b,paste0(d,"_minus_baseline"))
    if(loc=="kernel_gaussian") compare(a,b,paste0(d,"_minus_baseline"),"direct_squared")
  }
  if(d=="C1") compare(a,pick(design="C0"),"C1_minus_C0")
  if(d=="C2" && o==1) compare(a,pick(order=.4),"order_one_minus_point_four")
  if(d=="C0" && n==531441) compare(a,pick(size=65536),"M3_endpoint_minus_M2_endpoint")
}
bind <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
tables <- list(raw_metrics=bind(metrics),qq=bind(qq),diagnostic_moments=bind(dgstats),
  gram_conditional=bind(gram),paired_contrasts=bind(pairs))
projection <- list()
for(M in 1:3) for(design in c("C0","C1","C2","C0H","C2S")) {
  sch <- affine_schedule(if(M==1) 2048 else if(M==2) 4096 else 531441,design_id=design)
  fit <- affine_population_fit(sch,design);fresh <- affine_fresh(fit,c(.5,.5),design)
  projection[[length(projection)+1L]] <- data.frame(design_id=design,M=M,K=sch$K,
    sine_residual_variance=mean(affine_sine_cells(M)$variance),
    propensity_MSE=fresh$mse["pi"],cate_MSE_sum=sum(fresh$mse[c("pi","mu0","mu1")]),
    slope_MSE_sum=sum(fresh$mse[c("pi","m")]),four_MSE_sum=sum(fresh$mse),
    cate_population_signal_bias=fresh$cate$signal_bias,
    slope_population_root_bias=-fresh$slope$score_mean/(31/128+fresh$slope$coefficient_increment),
    reference_kind="analytic unpenalized projection, not a fitted estimator or contract supremum")
}
tables$projection_references <- bind(projection)
for(k in names(tables)) write.csv(tables[[k]],file.path(out,paste0(k,".csv")),row.names=FALSE)
manifest <- list(status="analysis_of_audited_saved_production",estimation_performed=FALSE,
  source_audit_sha256=digest::digest(file=file.path(accepted,"audit.json"),algo="sha256"),
  extractor_sha256=digest::digest(file=script,algo="sha256"),native_inputs=receipts,
  quantile_type=8,empirical_center="group mean of successful estimates, diagnostic only",
  paired_sign="left minus right",duplicate_policy="retain dnn profile; retain unique C0 n8192 dnn rows from kernel profile",
  data_hashes=lapply(names(tables),function(k) list(file=paste0(k,".csv"),rows=nrow(tables[[k]]),
    sha256=digest::digest(file=file.path(out,paste0(k,".csv")),algo="sha256"))))
jsonlite::write_json(manifest,file.path(out,"analysis_manifest.json"),pretty=TRUE,auto_unbox=TRUE,digits=17,na="null")
cat("Saved production analysis tables; no simulation or estimation performed.\n")
