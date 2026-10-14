# One estimator row per training/localization pair; SE routes share that point.
affine_quantile <- function(x,p) {
  x <- x[is.finite(x)]; if(!length(x)) return(rep(NA_real_,length(p)))
  unname(quantile(x,p,type=8))
}

affine_binomial_interval <- function(k,n) {
  if(n==0) return(c(NA_real_,NA_real_))
  c(if(k==0) 0 else qbeta(.025,k,n-k+1),if(k==n) 1 else qbeta(.975,k+1,n-k))
}

affine_coverage_comparison <- function(delta,length_gap) {
  n <- length(delta)
  if(length(length_gap)!=n || any(!delta%in%c(-1,0,1)) || any(!is.finite(length_gap)))
    stop("Invalid paired interval differences")
  data.frame(coverage_difference=if(n) mean(delta) else NA_real_,
    coverage_difference_mcse=if(n>1) sd(delta)/sqrt(n) else NA_real_,
    length_difference=if(n) mean(length_gap) else NA_real_,
    length_difference_mcse=if(n>1) sd(length_gap)/sqrt(n) else NA_real_,
    coverage_gains=sum(delta==1),coverage_losses=sum(delta== -1),
    coverage_discordances=sum(delta!=0),
    zero_discordance_one_sided_upper=if(n && all(delta==0)) 1-.05^(1/n) else NA_real_)
}

affine_training_pairs <- function(raw,cells) {
  groups <- c("cell_id","application","localization_method"); rows <- list()
  for(ix in split(seq_len(nrow(raw)),interaction(raw[groups],drop=TRUE,lex.order=TRUE))) {
    z <- raw[ix,,drop=FALSE]; a <- z[z$training_method=="full",,drop=FALSE]
    b <- z[z$training_method=="pooled_five_fold",,drop=FALSE]
    if(!nrow(a)) next
    b <- b[match(a$replication_id,b$replication_id),,drop=FALSE]
    planned <- cells$replications[match(a$cell_id[1],cells$cell_id)]
    point <- a$status=="ok" & b$status=="ok" & is.finite(a$estimate) & is.finite(b$estimate)
    point[is.na(point)] <- FALSE
    gap <- a$estimate[point]-b$estimate[point]
    loss <- (a$estimate[point]-a$point_target[point])^2-(b$estimate[point]-b$point_target[point])^2
    for(route in c("jackknife","direct_moment","direct_squared")) {
      available <- function(d) d[[paste0("variance_status_",route)]]!="unavailable" &
        is.finite(d[[paste0("SE_",route)]]) & d[[paste0("SE_",route)]]>=0 &
        is.finite(d[[paste0("lower_",route)]]) & is.finite(d[[paste0("upper_",route)]])
      valid <- point & available(a) & available(b); valid[is.na(valid)] <- FALSE
      cover <- function(d) as.integer(d[[paste0("lower_",route)]][valid]<=d$point_target[valid] &
        d[[paste0("upper_",route)]][valid]>=d$point_target[valid])
      dl <- (a[[paste0("upper_",route)]][valid]-a[[paste0("lower_",route)]][valid])-
        (b[[paste0("upper_",route)]][valid]-b[[paste0("lower_",route)]][valid])
      rows[[length(rows)+1L]] <- cbind(a[1,groups,drop=FALSE],data.frame(
        comparison="full_minus_pooled_five_fold",variance_method=route,planned=planned,
        point_pairs=sum(point),excluded_points=planned-sum(point),
        mean_gap=if(length(gap)) mean(gap) else NA_real_,
        RMS_gap=if(length(gap)) sqrt(mean(gap^2)) else NA_real_,
        MSE_difference=if(length(loss)) mean(loss) else NA_real_,
        MSE_difference_mcse=if(length(loss)>1) sd(loss)/sqrt(length(loss)) else NA_real_,
        paired=sum(valid),excluded=planned-sum(valid)),affine_coverage_comparison(cover(a)-cover(b),dl))
    }
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}

affine_boot_metrics <- function(est,target,se,ids,planned,B,seed) {
  if(B<2L) return(c(sd_mcse=NA_real_,rmse_mcse=NA_real_,ratio_mcse=NA_real_,q025_mcse=NA_real_,q975_mcse=NA_real_))
  # Resetting this deterministic seed gives common ID resamples for every method.
  ix <- affine_with_seed(seed,matrix(sample.int(planned,planned*B,replace=TRUE),planned,B))
  calc <- function(ii) {
    positions <- match(ii,ids); e <- est[positions]; t <- target[positions]; s <- se[positions]
    good <- is.finite(e); e <- e[good]; t <- t[good]; s <- s[good]
    sdv <- if(length(e)>1) sd(e) else NA_real_
    stu <- (e-t)/s; stu <- stu[is.finite(stu) & s>0]
    c(sd=sdv,rmse=if(length(e)) sqrt(mean((e-t)^2)) else NA_real_,
      ratio=if(is.finite(sdv)&&sdv>0&&any(is.finite(s))) sqrt(mean(s[is.finite(s)]^2))/sdv else NA_real_,
      affine_quantile(stu,c(.025,.975)))
  }
  draws <- vapply(seq_len(B),function(b) calc(ix[,b]),numeric(5))
  out <- apply(draws,1,function(v) if(sum(is.finite(v))>1) sd(v[is.finite(v)]) else NA_real_)
  setNames(out,c("sd_mcse","rmse_mcse","ratio_mcse","q025_mcse","q975_mcse"))
}

affine_boot_gap <- function(gap,ids,planned,B,seed) {
  if(B<2 || !length(gap)) return(c(sd_mcse=NA_real_,q95_absolute_mcse=NA_real_,median_mcse=NA_real_))
  ix <- affine_with_seed(seed,matrix(sample.int(planned,planned*B,replace=TRUE),planned,B))
  samples <- vapply(seq_len(B),function(b) {
    g <- gap[match(ix[,b],ids)]; g <- g[is.finite(g)]
    c(if(length(g)>1) sd(g) else NA_real_,affine_quantile(abs(g),.95),affine_quantile(g,.5))
  },numeric(3))
  setNames(apply(samples,1,function(g) if(sum(is.finite(g))>1) sd(g[is.finite(g)]) else NA_real_),
       c("sd_mcse","q95_absolute_mcse","median_mcse"))
}

# Across-localizer comparisons answer a different question from within-point
# variance comparisons. Positive differences are kernel minus DNN; all terms
# in each reported difference use the same replication intersection.
affine_localizer_pairs <- function(raw,cells) {
  if(!all(c("dnn","kernel_gaussian")%in%raw$localization_method)) return(data.frame())
  groups <- c("cell_id","application","training_method")
  rows <- list()
  for(ix in split(seq_len(nrow(raw)),interaction(raw[groups],drop=TRUE,lex.order=TRUE))) {
    z <- raw[ix,,drop=FALSE]; a <- z[z$localization_method=="dnn",,drop=FALSE]
    b <- z[z$localization_method=="kernel_gaussian",,drop=FALSE]
    if(!nrow(a)) next
    b <- b[match(a$replication_id,b$replication_id),,drop=FALSE]
    cell <- cells[match(a$cell_id[1],cells$cell_id),]
    planned <- if(a$training_method[1]%in%c("population_projection","population_ridge")) cell$control_replications else cell$replications
    point <- a$status=="ok" & b$status=="ok" & is.finite(a$estimate) & is.finite(b$estimate)
    point[is.na(point)] <- FALSE
    gap <- b$estimate[point]-a$estimate[point]
    loss <- (b$estimate[point]-b$point_target[point])^2-(a$estimate[point]-a$point_target[point])^2
    mean_or_na <- function(v) if(length(v)) mean(v) else NA_real_
    mcse <- function(v) if(length(v)>1) sd(v)/sqrt(length(v)) else NA_real_
    for(route in c("jackknife","direct_moment","direct_squared")) {
      available <- function(d) d[[paste0("variance_status_",route)]]!="unavailable" &
        is.finite(d[[paste0("SE_",route)]]) & d[[paste0("SE_",route)]]>=0 &
        is.finite(d[[paste0("lower_",route)]]) & is.finite(d[[paste0("upper_",route)]])
      valid <- point & available(a) & available(b); valid[is.na(valid)] <- FALSE
      coverage <- function(d) d[[paste0("lower_",route)]][valid]<=d$point_target[valid] &
        d[[paste0("upper_",route)]][valid]>=d$point_target[valid]
      delta <- as.integer(coverage(b))-as.integer(coverage(a))
      length_gap <- (b[[paste0("upper_",route)]][valid]-b[[paste0("lower_",route)]][valid])-
        (a[[paste0("upper_",route)]][valid]-a[[paste0("lower_",route)]][valid])
      rows[[length(rows)+1L]] <- cbind(a[1,groups,drop=FALSE],data.frame(
        comparison="kernel_minus_dnn",variance_method=route,planned=planned,
        point_pairs=sum(point),excluded_points=planned-sum(point),
        mean_estimate_difference=mean_or_na(gap),mean_estimate_difference_mcse=mcse(gap),
        mean_scaled_difference=mean_or_na(sqrt(a$N_reference[point])*gap),
        MSE_difference=mean_or_na(loss),MSE_difference_mcse=mcse(loss),
        mean_N_dnn=mean_or_na(a$N[point]),mean_N_kernel=mean_or_na(b$N[point]),
        mean_radius_dnn=mean_or_na(a$weighted_radius[point]),mean_radius_kernel=mean_or_na(b$weighted_radius[point]),
        interval_pairs=sum(valid),excluded_intervals=planned-sum(valid)),affine_coverage_comparison(delta,length_gap))
    }
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}

affine_diagnostic_row <- function(row,d) {
  scalar <- function(v) if(is.null(v)||!length(v)) NA_real_ else unname(v[1])
  out <- row[,c("cell_id","replication_id","application","training_method","localization_method","status"),drop=FALSE]
  accuracy <- d$accuracy
  if(!is.null(accuracy$fresh_functions)) {
    fresh <- accuracy$fresh_functions
    f <- fresh[[1]]
    f$mse <- accuracy$fold_average
    f$component_mse <- Reduce(`+`,lapply(fresh,`[[`,"component_mse"))/5
    for(app in c("cate","slope")) for(key in names(f[[app]]))
      f[[app]][[key]] <- mean(vapply(fresh,function(z) z[[app]][[key]],numeric(1)))
    f$integration_error <- max(vapply(fresh,`[[`,numeric(1),"integration_error"))
  } else f <- accuracy
  for(j in c("pi","mu0","mu1","m")) {
    out[[paste0("fresh_",j)]] <- scalar(f$mse[j])
    out[[paste0("fresh_local_mass_",j)]] <- scalar(accuracy$local_mass_average[j])
    out[[paste0("training_",j)]] <- scalar(d$training_mse[j])
    out[[paste0("clip_",j)]] <- scalar(d$clip_frequency[j])
  }
  app <- row$application
  out$mean_score_at_point <- scalar(f[[app]]$score_mean)
  out$coefficient_increment <- if(app=="slope") scalar(f$slope$coefficient_increment) else 0
  out$nuisance_product <- if(app=="slope") scalar(f$slope$product) else NA_real_
  out$local_target_mean_remainder <- if(app=="slope") scalar(f$slope$local_remainder) else NA_real_
  out$fixed_target_mixed_mean <- if(app=="slope") scalar(f$slope$mixed_spatial) else NA_real_
  out$integration_error <- scalar(f$integration_error)
  for(key in c("cate_b","slope_B","slope_b","slope_score")) out[[paste0("component_mse_",key)]] <- scalar(f$component_mse[key])
  for(key in c("transfer","component","target","cross","point_oracle_scaled","oracle_root_scaled","oracle_Gamma"))
    out[[paste0("covariance_",key)]] <- scalar(d$covariance[[app]][[key]])
  for(key in c("scaled_mean","scaled_centered","coefficient_spatial","identity_error",
      "weighted_coefficient_increment","weighted_product","weighted_local_remainder","weighted_mixed_spatial"))
    out[[paste0("decomposition_",key)]] <- scalar(d$decomposition[[app]][[key]])
  for(key in c("weighted_radius","radius90","radius95","largest_weight","mathematical_support","numerical_support","omitted_mass_bound","normalization_error"))
    out[[key]] <- scalar(d$local[[key]])
  out$grid_prediction_radius_squared <- scalar(d$deletions$grid_prediction_radius_squared)
  out$deletion_order_reference_K_squared_over_n <- scalar(d$deletions$order_reference_K_squared_over_n)
  out$scaled_deletion_order_reference <- scalar(d$deletions$scaled_order_reference)
  out$population_spatial_drift <- scalar(d$population_localizer$drift)
  out$population_integration_error <- scalar(d$population_localizer$integration_error)
  out$conditional_noise_VT <- scalar(d$conditional_noise[paste0(app,"_VT")])
  out$kernel_deletion_recomputations <- scalar(d$kernel_deletion_recomputations)
  out$fitted_score_scaled <- row$score_cov_scaled_jackknife
  out$fitted_target_scaled <- row$N*row$variance_jackknife
  for(scale in c("score","target")) for(ref in c("point_oracle","oracle_root","asymptotic")) {
    name <- paste("ratio",scale,ref,sep="_")
    numerator <- if(scale=="score") out$fitted_score_scaled else out$fitted_target_scaled
    denominator <- if(ref=="asymptotic") {
      if(app=="cate") 1/(4*sqrt(15)) else if(scale=="score") 31/8192 else 2/31
    } else {
      value <- out[[paste0("covariance_",ref,"_scaled")]]
      if(scale=="target") value/out$covariance_oracle_Gamma^2 else value
    }
    out[[paste0(name,"_denominator")]] <- denominator
    out[[paste0(name,"_denominator_status")]] <- if(!is.finite(denominator)) "unavailable" else
      if(denominator==0) "zero" else if(denominator<0) "negative" else "positive"
    out[[name]] <- if(is.finite(numerator) && is.finite(denominator) && denominator>0) numerator/denominator else NA_real_
  }
  out
}

# Keep the single-row diagnostic definition, but release each bounded batch
# after binding it. A production run must not retain a million tiny frames.
affine_diagnostic_table <- function(raw,diagnostics,block_size=256L) {
  if(length(block_size)!=1L || !is.finite(block_size) || block_size<1 || block_size!=floor(block_size))
    stop("Invalid diagnostic batch size")
  if(!nrow(raw)) return(data.frame())
  batches <- split(seq_len(nrow(raw)),ceiling(seq_len(nrow(raw))/block_size))
  parts <- lapply(batches,function(indices) {
    rows <- lapply(indices,function(i) {
      key <- affine_diagnostic_key(raw$cell_id[i],raw$replication_id[i],raw$training_method[i],raw$localization_method[i])
      affine_diagnostic_row(raw[i,,drop=FALSE],diagnostics[[key]])
    })
    do.call(rbind,rows)
  })
  do.call(rbind,unname(parts))
}

affine_summarise <- function(raw,diagnostics,cell_table,bootstrap=199L,bootstrap_seed=20261015L) {
  if(!nrow(raw)) return(c(setNames(rep(list(data.frame()),10L),c("performance","paired","variance_pairs",
    "diagnostics","gram_events","spatial_decomposition","covariance_ratios","error_counts","localizer_pairs","training_pairs")),
    list(bootstrap_replications=bootstrap,bootstrap_seed=bootstrap_seed,quantile_type=8L)))
  group_keys <- c("cell_id","application","training_method","localization_method")
  keys <- interaction(raw[group_keys],drop=TRUE,lex.order=TRUE)
  routes <- c("jackknife","direct_moment","direct_squared")
  performance <- list(); variance_pairs <- list(); paired <- list()
  for(ix in split(seq_len(nrow(raw)),keys)) {
    z <- raw[ix,,drop=FALSE]; cell <- cell_table[match(z$cell_id[1],cell_table$cell_id),]
    planned <- if(z$training_method[1]%in%c("population_projection","population_ridge")) cell$control_replications else cell$replications
    if(is.null(planned)) planned <- max(z$replication_id)
    good <- z$status=="ok" & is.finite(z$estimate)
    e <- z$estimate[good]; target <- z$point_target[good]
    mcSD <- if(length(e)>1) sd(e) else NA_real_
    idfields <- z[1,group_keys,drop=FALSE]
    cover <- avail <- lengths <- list()
    for(route in routes) {
      se <- z[[paste0("SE_",route)]]
      lower <- z[[paste0("lower_",route)]]; upper <- z[[paste0("upper_",route)]]
      valid <- good & z[[paste0("variance_status_",route)]]!="unavailable" & is.finite(se) & se>=0 & is.finite(lower) & is.finite(upper)
      covered <- valid & lower<=z$point_target & upper>=z$point_target
      cover[[route]] <- covered; avail[[route]] <- valid; lengths[[route]] <- upper-lower
      intervals <- sum(valid); cov <- if(intervals) sum(covered)/intervals else NA_real_
      ci <- affine_binomial_interval(sum(covered),intervals)
      stumask <- valid & se>0; student <- (z$estimate[stumask]-z$point_target[stumask])/se[stumask]
      qq <- affine_quantile(student,c(.025,.5,.975))
      meanse <- if(any(valid)) mean(se[valid]) else NA_real_
      calibration_SD <- if(sum(valid)>1) sd(z$estimate[valid]) else NA_real_
      vals <- list(variance_method=route,
         variance_role=if(z$localization_method[1]=="kernel_gaussian" && route=="direct_squared")
           "conventional_kernel_sandwich" else if(route=="jackknife") "frozen_score_jackknife" else "direct_diagnostic",
         planned=planned,completed=nrow(z),successful_estimates=sum(good),
         intervals=intervals,studentized=sum(stumask),failures=sum(!good),degenerate=sum(valid & se==0),
         bias=if(length(e)) mean(e-target) else NA_real_,sd=mcSD,rmse=if(length(e)) sqrt(mean((e-target)^2)) else NA_real_,
         bias_over_SD=if(length(e) && is.finite(mcSD) && mcSD>0) mean(e-target)/mcSD else NA_real_,
         calibration_pairs=sum(valid),calibration_SD=calibration_SD,
         mean_variance_over_MC_variance=if(is.finite(calibration_SD) && calibration_SD>0)
           mean(se[valid]^2)/calibration_SD^2 else NA_real_,
         mean_N=if(any(good)) mean(z$N[good]) else NA_real_,N_reference=z$N_reference[1],
         mean_weighted_radius=if(any(good)) mean(z$weighted_radius[good]) else NA_real_,
         mean_largest_weight=if(any(good)) mean(z$largest_weight[good]) else NA_real_,
         bias_mcse=mcSD/sqrt(sum(good)),coverage=cov,coverage_mcse=if(intervals) sqrt(cov*(1-cov)/intervals) else NA_real_,
         coverage_lower=ci[1],coverage_upper=ci[2],planned_covering_fraction=sum(covered)/planned,
         lower_miss=if(intervals) sum(valid & upper<z$point_target)/intervals else NA_real_,
         upper_miss=if(intervals) sum(valid & lower>z$point_target)/intervals else NA_real_,
         mean_SE=meanse,mean_SE_over_SD=if(is.finite(mcSD)&&mcSD>0) meanse/mcSD else NA_real_,
         rms_SE_over_SD=if(is.finite(mcSD)&&mcSD>0&&intervals) sqrt(mean(se[valid]^2))/mcSD else NA_real_,
         mean_length=if(intervals) mean(lengths[[route]][valid]) else NA_real_,
         student_mean=if(length(student)) mean(student) else NA_real_,student_sd=if(length(student)>1) sd(student) else NA_real_,
         student_q025=qq[1],student_median=qq[2],student_q975=qq[3],
         student_lower_tail=if(length(student)) mean(student < -qnorm(.975)) else NA_real_,
         student_upper_tail=if(length(student)) mean(student > qnorm(.975)) else NA_real_)
      bs <- affine_boot_metrics(ifelse(good,z$estimate,NA_real_),z$point_target,
            ifelse(valid,se,NA_real_),z$replication_id,planned,bootstrap,bootstrap_seed)
      performance[[length(performance)+1L]] <- cbind(idfields,as.data.frame(c(vals,as.list(bs))))
    }
    for(comp in list(c("jackknife","direct_moment"),c("jackknife","direct_squared"))) {
      valid <- avail[[comp[1]]]&avail[[comp[2]]]
      delta <- as.integer(cover[[comp[1]]][valid])-as.integer(cover[[comp[2]]][valid])
      dl <- lengths[[comp[1]]][valid]-lengths[[comp[2]]][valid]
      variance_pairs[[length(variance_pairs)+1L]] <- cbind(idfields,data.frame(comparison=paste(comp,collapse="_vs_"),
          paired=sum(valid),excluded=planned-sum(valid)),affine_coverage_comparison(delta,dl))
    }
  }
  for(cell_id in unique(raw$cell_id)) for(app in c("cate","slope"))
      for(localization_method in unique(raw$localization_method)) for(method in c("pooled_five_fold","oracle")) {
    a <- raw[raw$cell_id==cell_id & raw$application==app & raw$training_method=="full" & raw$localization_method==localization_method,,drop=FALSE]
    b <- raw[raw$cell_id==cell_id & raw$application==app & raw$training_method==method & raw$localization_method==localization_method,,drop=FALSE]
    plan <- cell_table$replications[match(cell_id,cell_table$cell_id)]
    pos <- match(a$replication_id,b$replication_id)
    valid <- a$status=="ok" & is.finite(a$estimate) & !is.na(pos)
    valid <- valid & !is.na(b$status[pos]) & b$status[pos]=="ok" & is.finite(b$estimate[pos])
    valid[is.na(valid)] <- FALSE; ar <- a[valid,]; br <- b[pos[valid],]
    if(nrow(ar) && any(ar$N!=br$N)) stop("Paired methods have different effective sizes")
    gap <- sqrt(ar$N)*(ar$estimate-br$estimate)
    Omega <- if(app=="cate") 1/(4*sqrt(15)) else 2/31
    qq <- affine_quantile(abs(gap),.95)
    paired[[length(paired)+1L]] <- data.frame(cell_id=cell_id,application=app,localization_method=localization_method,
       comparison=paste0("full_vs_",if(method=="oracle") "oracle" else "CF"),
       planned=plan,paired=length(gap),excluded=plan-length(gap),mean=if(length(gap)) mean(gap) else NA_real_,
       mean_mcse=if(length(gap)>1) sd(gap)/sqrt(length(gap)) else NA_real_,sd=if(length(gap)>1) sd(gap) else NA_real_,
       median=affine_quantile(gap,.5),q95_absolute=qq,mean_Omega_units=if(length(gap)) mean(gap)/sqrt(Omega) else NA_real_,
       q95_absolute_Omega_units=qq/sqrt(Omega))
    bs <- affine_boot_gap(gap,ar$replication_id,plan,bootstrap,bootstrap_seed)
    paired[[length(paired)]] <- cbind(paired[[length(paired)]],as.data.frame(as.list(bs)))
  }
  # Raw diagnostic rows preserve receipt keys; scientific interpretation remains separate.
  dg <- affine_diagnostic_table(raw,diagnostics)
  events <- list()
  for(cell_id in unique(raw$cell_id)) {
    z <- raw[raw$cell_id==cell_id & raw$application=="cate" & raw$training_method=="full" & raw$localization_method=="dnn",,drop=FALSE]
    observed <- !is.na(z$gram_event); R <- sum(observed); k <- sum(!z$gram_event[observed])
    ci <- affine_binomial_interval(k,R)
    events[[length(events)+1L]] <- data.frame(cell_id=cell_id,completed=nrow(z),observed=R,unknown=sum(!observed),
      failure_fraction=if(R) k/R else NA_real_,scaled_failure=if(R) z$N[1]*k/R else NA_real_,
      lower=ci[1],upper=ci[2],zero_failure_one_sided_upper=if(R && k==0) 1-.05^(1/R) else NA_real_)
  }
  bind <- function(rows) if(length(rows)) do.call(rbind,rows) else data.frame()
  c(list(performance=bind(performance),paired=bind(paired),variance_pairs=bind(variance_pairs),
       diagnostics=dg,gram_events=bind(events),localizer_pairs=affine_localizer_pairs(raw,cell_table),
       training_pairs=affine_training_pairs(raw,cell_table)),
       affine_extra_summaries(raw,dg,cell_table),
       list(bootstrap_replications=bootstrap,bootstrap_seed=bootstrap_seed,quantile_type=8L))
}


# These summaries use the same successful replication intersection for every
# term in a decomposition, rather than adding variances with different masks.
affine_extra_summaries <- function(raw,diagnostic_rows,cells) {
  groups <- c("cell_id","application","training_method","localization_method")
  split_rows <- split(seq_len(nrow(raw)),interaction(raw[groups],drop=TRUE,lex.order=TRUE))
  spatial <- ratios <- errors <- list()
  safe_mean <- function(v) if(length(v)) mean(v) else NA_real_
  safe_var <- function(v) if(length(v)>1) var(v) else NA_real_
  for(ix in split_rows) {
    z <- raw[ix,,drop=FALSE]; d <- diagnostic_rows[ix,,drop=FALSE]; ids <- z[1,groups,drop=FALSE]
    cell <- cells[match(z$cell_id[1],cells$cell_id),]
    planned <- if(z$training_method[1]%in%c("population_projection","population_ridge")) cell$control_replications else cell$replications
    good <- z$status=="ok" & is.finite(z$estimate) & is.finite(z$spatial_drift) & is.finite(z$centered_remainder)
    drift <- z$spatial_drift[good]; remainder <- z$centered_remainder[good]; total <- z$estimate[good]-z$point_target[good]
    covariance <- if(length(drift)>1) cov(drift,remainder) else NA_real_
    Omega <- if(z$application[1]=="cate") 1/(4*sqrt(15)) else 2/31
    spatial[[length(spatial)+1L]] <- cbind(ids,data.frame(planned=planned,paired=sum(good),excluded=planned-sum(good),
      mean_drift=safe_mean(drift),mean_remainder=safe_mean(remainder),
      drift_variance=safe_var(drift),remainder_variance=safe_var(remainder),drift_remainder_covariance=covariance,
      total_variance=safe_var(total),variance_identity_error=safe_var(total)-safe_var(drift)-safe_var(remainder)-2*covariance,
      mean_drift_Omega_units=safe_mean(sqrt(z$N[good])*drift/sqrt(Omega)),
      mean_drift_mcse=if(length(drift)>1) sd(drift)/sqrt(length(drift)) else NA_real_))
    for(scale in c("score","target")) for(ref in c("point_oracle","oracle_root","asymptotic")) {
      name <- paste("ratio",scale,ref,sep="_"); value <- d[[name]]
      denominator <- d[[paste0(name,"_denominator")]]; ds <- d[[paste0(name,"_denominator_status")]]
      numerator <- if(scale=="score") d$fitted_score_scaled else d$fitted_target_scaled
      valid <- z$status=="ok" & is.finite(value); count <- sum(valid)
      ratios[[length(ratios)+1L]] <- cbind(ids,data.frame(scale=scale,reference=ref,planned=planned,
        completed=nrow(z),paired=count,excluded=planned-count,positive_denominators=sum(ds=="positive"),
        zero_denominators=sum(ds=="zero"),negative_denominators=sum(ds=="negative"),unavailable_denominators=sum(ds=="unavailable"),
        unavailable_numerators=sum(!is.finite(numerator)),mean_ratio=safe_mean(value[valid]),
        mean_ratio_mcse=if(count>1) sd(value[valid])/sqrt(count) else NA_real_,
        median_ratio=affine_quantile(value[valid],.5),q95_ratio=affine_quantile(value[valid],.95),
        ratio_of_means=if(count) mean(numerator[valid])/mean(denominator[valid]) else NA_real_))
    }
    add_errors <- function(stage,route,status,category) {
      category[is.na(category)] <- "none"
      tab <- as.data.frame(table(status=status,category=category),stringsAsFactors=FALSE)
      tab <- tab[tab$Freq>0,,drop=FALSE]
      cbind(ids[rep(1,nrow(tab)),,drop=FALSE],data.frame(stage=stage,variance_method=route,
           status=as.character(tab$status),error_category=as.character(tab$category),count=tab$Freq))
    }
    errors[[length(errors)+1L]] <- add_errors("point","none",z$status,z$error_category)
    errors[[length(errors)+1L]] <- add_errors("diagnostic","none",z$diagnostic_status,
       ifelse(z$diagnostic_status=="failed","diagnostic_error",NA_character_))
    for(route in c("jackknife","direct_moment","direct_squared")) errors[[length(errors)+1L]] <-
      add_errors("variance",route,z[[paste0("variance_status_",route)]],
          ifelse(!is.na(z[[paste0("variance_error_",route)]]),"covariance_error",NA_character_))
  }
  list(spatial_decomposition=do.call(rbind,spatial),covariance_ratios=do.call(rbind,ratios),error_counts=do.call(rbind,errors))
}
