# New orders from validated saved signals; independent of nuisance fitting.
.sipp_sensitivity_sources <- Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile))
.sipp_sensitivity_file <- tail(.sipp_sensitivity_sources,1L)[[1L]]
sipp_401k_localization_sensitivity <- function(fit, orders, profile_ids = NULL) {
  n <- nrow(fit$analysis_data)
  if (!length(orders) || any(!is.finite(orders)) || any(orders < 1 | orders >= n) ||
      any(orders != floor(orders)) || anyDuplicated(orders))
    stop("Orders must be distinct integers strictly between zero and sample size.",call.=FALSE)
  baseline <- analyze_sipp_401k_frozen(fit)$information
  if (is.null(profile_ids)) profile_ids <- fit$profiles$profile_id
  if (!length(profile_ids) || anyDuplicated(profile_ids) ||
      any(!profile_ids %in% fit$profiles$profile_id)) stop("Invalid sensitivity profiles.",call.=FALSE)
  d <- fit$analysis_data; cols <- fit$specification$metric$columns
  if (is.null(cols)) cols <- c("age_z","asinh_inc_z")
  X <- as.matrix(d[,cols,drop=FALSE])
  output <- list()
  for (scale in names(fit$dnn_weights)) for (id in profile_ids) {
    j <- match(id,fit$profiles$profile_id)
    x0 <- as.numeric(fit$profiles[j,cols,drop=FALSE])
    distance <- sqrt(rowSums(sweep(X,2L,x0,"-")^2)); ord <- order(distance,seq_len(n))
    payload <- fit$dnn_weights[[scale]]
    mechanisms <- names(payload)[!startsWith(names(payload),"pseudo_outcome_")]
    for (mechanism in mechanisms) for (s in sort(orders)) {
      method <- paste0("dnn_",mechanism); phi <- payload[[paste0("pseudo_outcome_",mechanism)]]
      e <- fit$nuisances[[mechanism]]$e_hat
      z <- sipp_401k_frozen_rank(phi[ord],s); w <- numeric(n); w[ord] <- z$weights
      pc <- z$centered[,1L]; square <- sum(pc^2); h <- if(square>0)pc^2/square else rep(NA_real_,n)
      top <- if(square>0)which.max(h) else 1L
      base <- baseline[baseline$outcome_scale==scale & baseline$profile_id==id &
        baseline$method==method,]
      tracked <- match(base$top_row_id,d$row_id[ord])
      radius_index <- which(cumsum(z$weights)>=.9)[1L]; radius <- distance[ord][radius_index]
      in_radius <- distance <= radius
      output[[length(output)+1L]] <- data.frame(outcome_scale=scale,profile_id=id,
        method=method,s=s,tau=z$tau,se=sqrt(z$covariance[1L,1L]),ess=1/sum(w^2),
        radius_90_z=radius,treated_mass=sum(w*d$D),
        ipw_treated_kish=sipp_401k_kish_effective_size(w*d$D/e),
        ipw_control_kish=sipp_401k_kish_effective_size(w*(1-d$D)/(1-e)),
        jk_top1_sq_share=if(square>0)max(h) else NA_real_,
        jk_top5_sq_share=if(square>0)sum(head(sort(h,decreasing=TRUE),5L)) else NA_real_,
        top_row_id=d$row_id[ord[top]],frozen_delete_tau=z$loo[top,1L],
        nearest_weight=z$weights[1L],baseline_top_row_id=base$top_row_id,
        baseline_case_weight=z$weights[tracked],
        baseline_case_jk_sq_share=h[tracked],
        age_min_90=if("age"%in%names(d))min(d$age[in_radius]) else NA_real_,
        age_max_90=if("age"%in%names(d))max(d$age[in_radius]) else NA_real_,
        income_min_90=if("inc"%in%names(d))min(d$inc[in_radius]) else NA_real_,
        income_max_90=if("inc"%in%names(d))max(d$inc[in_radius]) else NA_real_)
    }
  }
  do.call(rbind,output)
}

export_sipp_401k_localization <- function(fit,orders,output_dir,allow_validation=FALSE) {
  if (dir.exists(output_dir)||file.exists(output_dir))
    stop("Sensitivity directory already exists.",call.=FALSE)
  validation <- isTRUE(fit$run_info$smoke) || (!is.null(fit$run_info$sample_n) &&
    !is.na(fit$run_info$sample_n)) || fit$source$n_fit != fit$source$n_analysis
  if (validation && !allow_validation)
    stop("Smoke/subsample sensitivity needs an explicit validation override.",call.=FALSE)
  if (!length(orders) || any(!is.finite(orders)) || any(orders < 1) ||
      any(orders >= nrow(fit$analysis_data)) || any(orders != floor(orders)) ||
      anyDuplicated(orders))
    stop("Orders must be distinct feasible integers.",call.=FALSE)
  orders <- sort(unique(c(orders,fit$estimates$s[is.finite(fit$estimates$s)])))
  rows <- sipp_401k_localization_sensitivity(fit,orders)
  if (!dir.create(output_dir,recursive=TRUE)) stop("Cannot create sensitivity directory.",call.=FALSE)
  write.csv(rows,file.path(output_dir,"localization_all.csv"),row.names=FALSE)
  selected <- c("age40_inc_p25","age40_inc_p50","age40_inc_p75")
  if (!all(selected%in%fit$profiles$profile_id)) selected <- head(fit$profiles$profile_id,3L)
  method <- if("dnn_oob"%in%fit$estimates$method)"dnn_oob" else "dnn_fullsample"
  for (scale in unique(rows$outcome_scale)) {
    r <- rows[rows$outcome_scale==scale & rows$method==method & rows$profile_id%in%selected,]
    multiplier <- switch(scale,raw_dollars=.001,capped_dollars=.001,positive_assets=100,asinh_dollars=1)
    units <- switch(scale,raw_dollars="$1,000",capped_dollars="capped $1,000",
                    positive_assets="percentage points",asinh_dollars="asinh dollars")
    z <- qnorm(1-fit$config$alpha/2)
    limits <- range(c(r$tau-z*r$se,r$tau+z*r$se))*multiplier
    if(diff(limits)==0)limits<-c(-1,1)
    pdf(file.path(output_dir,paste0(scale,"_localization.pdf")),width=10.5,height=4.4)
    tryCatch({
      par(mfrow=c(1,length(selected)),mar=c(5,4,3,1),oma=c(2,0,0,0),las=1)
      for(id in selected) {
        a <- r[r$profile_id==id,]; income <- fit$profiles$inc[match(id,fit$profiles$profile_id)]
        plot(a$s,a$tau*multiplier,type="b",pch=16,col="#17688c",ylim=limits,
          xlab="Exact DNN order s",ylab=paste0("Contrast / interval (",units,")"),
          main=paste0(if(validation)"VALIDATION: " else "", "Income $",
                     format(round(income),big.mark=",",scientific=FALSE,trim=TRUE)),bty="l")
        arrows(a$s,(a$tau-z*a$se)*multiplier,a$s,(a$tau+z*a$se)*multiplier,
          angle=90,code=3,length=.035,col="#17688c")
        abline(h=0,col="grey60")
        abline(v=unique(fit$estimates$s),lty=3,col="#a65d28")
      }
      mtext("Saved single-fit/OOB signals held fixed; dotted line marks the retained default.\nSmaller s broadens support. This comparison does not estimate spatial bias.",
        side=1,outer=TRUE,line=.2,cex=.8)
    },finally=dev.off())
  }
  files <- list.files(output_dir,full.names=TRUE);names(files)<-basename(files)
  metadata <- list(purpose=if(validation)"validation" else "empirical",
    fit_fingerprint=sipp_401k_fit_fingerprint(fit),signals_frozen=TRUE,
    orders=orders,retained_default=unique(fit$estimates$s),all_profiles=fit$profiles$profile_id,
    displayed_profiles=selected,displayed_method=method,
    interpretation="Frozen localization sensitivity; no new nuisance fit or selected replacement default.",
    analysis_source_sha256=sipp_401k_file_hash(.sipp_sensitivity_file),
    files=lapply(files,sipp_401k_file_hash))
  jsonlite::write_json(metadata,file.path(output_dir,"sensitivity.json"),pretty=TRUE,
    auto_unbox=TRUE,null="null",na="null",digits=NA)
  invisible(metadata)
}
