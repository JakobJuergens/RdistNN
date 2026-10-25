# Manuscript publication from pinned portable SIPP evidence; no learner fitting.
.sipp_manuscript_frames <- Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile))
.sipp_manuscript_source <- normalizePath(tail(.sipp_manuscript_frames,1L)[[1L]])
.sipp_localization_source <- file.path(dirname(.sipp_manuscript_source),"localization_sensitivity.R")
.sipp_affine_source <- file.path(dirname(.sipp_manuscript_source),"affine_export.R")
source(.sipp_localization_source)
source(.sipp_affine_source)

render_sipp_401k_manuscript <- function(input_dir,historical_dir,output_dir,table_dir,
                                      repo_root,check_only=FALSE) {
  for(pkg in c("digest","jsonlite"))
    if(!requireNamespace(pkg,quietly=TRUE))stop("Required package: ",pkg,call.=FALSE)
  if(!check_only&&!capabilities("cairo"))
    stop("SIPP PDF rendering requires an R installation with Cairo support.",call.=FALSE)
  sha <- function(path)digest::digest(file=path,algo="sha256")
  json <- function(path)jsonlite::read_json(path,simplifyVector=TRUE)
  repo_root <- normalizePath(repo_root,winslash="/")
  relative <- function(path) {
    p <- normalizePath(path,winslash="/",mustWork=FALSE)
    if(startsWith(p,paste0(repo_root,"/")))substring(p,nchar(repo_root)+2L) else p
  }
  inputs <- list()
  consume <- function(path,expected=NULL) {
    if(!file.exists(path))stop("Required input missing: ",path,call.=FALSE)
    h <- sha(path)
    if(!is.null(expected)&&!identical(h,as.character(expected)))
      stop("Input hash mismatch: ",relative(path),call.=FALSE)
    inputs[[relative(path)]] <<- h
    path
  }
  receipt <- json(consume(file.path(input_dir,"receipt.json")))
  for(name in names(receipt$files))consume(file.path(input_dir,name),receipt$files[[name]])
  evidence <- file.path(input_dir,"evidence")
  selection <- json(file.path(evidence,"selection.json"))
  sensitivity <- json(file.path(evidence,"sensitivity.json"))
  accepted <- json(file.path(evidence,"acceptance_verification.json"))
  fit <- json(file.path(evidence,"fit_manifest.json"))
  averaged_ties <- identical(accepted$localization_policy,"exact_distance_block_averaging_v1")
  if(averaged_ties)stopifnot(isTRUE(accepted$averaged_weights_reconstructed),
    isTRUE(accepted$deletions_reconstructed),isTRUE(accepted$represented_distances_verified))
  if(selection$purpose!="empirical"||sensitivity$purpose!="empirical")
    stop("Manuscript publication requires empirical evidence; validation is refused.",call.=FALSE)
  stopifnot(identical(selection$fit_fingerprint,sensitivity$fit_fingerprint),
    identical(selection$fit_fingerprint,accepted$fit_fingerprint),
    accepted$n==9915L,fit$status=="complete",sensitivity$signals_frozen,
    sensitivity$retained_default==397L,accepted$old_oob_signal_identical,
    accepted$old_crossfit_signal_identical,accepted$common_final_weights_identical,
    fit$specification$adjustment=="baseline",fit$specification$reporting=="age_income")
  historic <- json(consume(file.path(historical_dir,"provenance.json")))
  historic_path <- consume(file.path(historical_dir,"accepted_dnn_estimates.csv"),
    historic$evidence_sha256[["accepted_dnn_estimates.csv"]])
  read <- function(name)read.csv(file.path(evidence,paste0(name,".csv")),stringsAsFactors=FALSE)
  estimates <- read("estimates_all"); info <- read("information"); support <- read("support")
  paired <- read("paired_training"); means <- read("component_means")
  covariance <- read("component_covariance"); orders <- read("localization_all")
  old <- read.csv(historic_path,stringsAsFactors=FALSE)
  methods <- c("dnn_fullsample","dnn_oob","dnn_crossfit")
  profiles <- paste0("age40_inc_p",c("05","10","25","50","75","90","95"))
  display <- profiles[3:5]
  grid <- c(99L,149L,198L,248L,298L,397L,496L,745L)
  keys <- function(x,fields)do.call(paste,c(x[fields],sep="|"))
  stopifnot(nrow(estimates)==21L,!anyDuplicated(keys(estimates,c("profile_id","method"))),
    setequal(estimates$method,methods),setequal(estimates$profile_id,profiles),
    all(estimates$outcome_scale=="raw_dollars"),all(estimates$s==397L),
    setequal(selection$selected_profiles,display),setequal(sensitivity$orders,grid),
    nrow(orders)==168L,!anyDuplicated(keys(orders,c("profile_id","method","s"))),
    setequal(orders$s,grid),setequal(orders$method,methods),setequal(orders$profile_id,profiles),
    all(orders$outcome_scale=="raw_dollars"),nrow(means)==63L,nrow(covariance)==189L,
    nrow(paired)==21L,nrow(info)==21L,nrow(support)==21L)
  for(m in methods)stopifnot(setequal(estimates$profile_id[estimates$method==m],profiles))
  one <- function(x,id,method=NULL,s=NULL) {
    keep <- x$profile_id==id
    if(!is.null(method))keep <- keep&x$method==method
    if(!is.null(s))keep <- keep&x$s==s
    r <- x[keep,,drop=FALSE]
    if(nrow(r)!=1L)stop("Incomplete or duplicate evidence selection.",call.=FALSE)
    r
  }
  close <- function(a,b)stopifnot(isTRUE(all.equal(a,b,tolerance=1e-8,check.attributes=FALSE)))
  for(id in profiles)for(m in methods) {
    e <- one(estimates,id,m); z <- one(info,id,m)
    c <- covariance[covariance$profile_id==id&covariance$method==m,,drop=FALSE]
    u <- means[means$profile_id==id&means$method==m,,drop=FALSE]
    stopifnot(nrow(c)==9L,nrow(u)==3L)
    close(sum(c$covariance),e$se^2);close(sum(u$mean),e$tau)
    close(z$tau,e$tau);close(z$se,e$se)
    for(s in grid)one(orders,id,m,s)
    if(m!="dnn_fullsample"&&!averaged_ties) {
      old_method <- if(m=="dnn_oob")"dnn_fullsample" else m
      h <- one(old[old$outcome_scale=="raw_dollars",],id,old_method)
      close(e$tau,h$tau);close(e$se,h$se)
    }
  }
  asinh <- old[old$outcome_scale=="asinh_dollars",]
  stopifnot(nrow(asinh)==14L,setequal(asinh$method,c("dnn_fullsample","dnn_crossfit")))
  labels <- c(dnn_fullsample="All-tree",dnn_oob="Single-fit/OOB",dnn_crossfit="Five-fold")
  plabel <- function(id)paste0(as.integer(sub("age40_inc_p","",id)),"\\%")
  value <- function(x,d=2L)formatC(round(x,d),format="f",digits=d,big.mark=",",decimal.mark=".")
  number <- function(x,d=2L)paste0("\\(",value(x,d),"\\)")
  interval <- function(lo,hi)paste0("\\([",value(lo),",\\ ",value(hi),"]\\)")
  estimate <- function(tau,se)paste0("\\(",value(tau/1000),"\\;(",value(se/1000),")\\)")
  line <- function(cells)paste0(paste(cells,collapse=" & ")," \\\\")
  tables <- list(); table_records <- list()
  add <- function(name,header,rows,input_rows=length(rows),columns=NULL) {
    if(is.null(columns))columns <- paste0("@{}",paste0("l",paste(rep("r",length(header)-1L),collapse="")),"@{}")
    tables[[paste0(name,".tex")]] <<- c("% Generated by 07_render_affine_manuscript.R; do not edit.",
      paste0("\\begin{tabular}{",columns,"}"),"\\hline",line(header),"\\hline",rows,
      "\\hline","\\end{tabular}")
    table_records[[name]] <<- list(rows=length(rows),input_rows=as.integer(input_rows),
      file=relative(file.path(table_dir,paste0(name,".tex"))))
  }
  rows <- list()
  for(id in profiles)for(m in methods) {
    e <- one(estimates,id,m); z <- one(info,id,m)
    rows[[length(rows)+1L]] <- line(c(plabel(id),labels[m],number(e$tau/1000),number(e$se/1000),
      interval(e$ci_lo/1000,e$ci_hi/1000),number(100*z$jk_top1_sq_share,1)))
  }
  add("estimates_all",c("Income","Prediction","Estimate","SE","95\\% interval","Largest (\\%)"),
      unlist(rows),21L,columns="@{}llrrrr@{}")
  rows <- vapply(profiles,function(id) {
    z <- one(info,id,"dnn_oob"); w <- one(support,id,"dnn_oob")
    line(c(plabel(id),number(w$radius_90_z,4),number(w$treated_mass,3),
      number(z$ipw_treated_kish,1),number(z$ipw_control_kish,1),number(100*z$jk_top5_sq_share,1)))
  },character(1))
  add("information",c("Income","Radius","Eligible mass","IPW eligible","IPW control","Top five (\\%)"),rows)
  rows <- vapply(profiles,function(id) {
    z <- one(info,id,"dnn_oob")
    line(c(plabel(id),number(z$top_D,0),number(z$top_Y/1000),number(z$top_propensity,3),
      number(z$top_rank,0),number(100*z$top_weight,2),number(z$frozen_delete_tau/1000)))
  },character(1))
  add("leading_cases",c("Income","Eligible","Assets","Propensity","Rank","Weight (\\%)","Deleted estimate"),rows)
  rows <- list()
  for(id in profiles)for(i in which(paired$profile_id==id)) {
    p <- paired[i,]; a <- one(estimates,id,p$method_1); b <- one(estimates,id,p$method_2)
    close(p$gap,b$tau-a$tau);close(p$paired_gap_se^2,a$se^2+b$se^2-2*p$covariance)
    rows[[length(rows)+1L]] <- line(c(plabel(id),paste0(labels[p$method_2]," $-$ ",labels[p$method_1]),
      number(p$gap/1000),number(p$paired_gap_se/1000),number(p$absolute_signal_movement/1000),
      number(p$correlation,4)))
  }
  add("paired",c("Income","Prediction difference","Gap","Gap SE","Absolute movement","Correlation"),
    unlist(rows),21L,columns="@{}llrrrr@{}")
  rows <- covrows <- list()
  components <- c("contrast","eligible_residual","ineligible_residual")
  for(id in profiles)for(m in methods) {
    a <- means[means$profile_id==id&means$method==m,]
    a <- a$mean[match(components,a$component)]
    rows[[length(rows)+1L]] <- line(c(plabel(id),labels[m],vapply(a/1000,number,character(1)),number(sum(a)/1000)))
    c <- covariance[covariance$profile_id==id&covariance$method==m,]
    at <- function(i,j)c$covariance[c$component_1==components[i]&c$component_2==components[j]]
    v <- c(at(1,1),at(2,2),at(3,3),at(1,2)+at(2,1),at(1,3)+at(3,1),at(2,3)+at(3,2),sum(c$covariance))
    covrows[[length(covrows)+1L]] <- line(c(plabel(id),labels[m],vapply(v/1e6,function(x)number(x,1),character(1))))
  }
  add("component_means",c("Income","Prediction","Contrast","Eligible correction","Control correction","Total"),
      unlist(rows),63L,columns="@{}llrrrr@{}")
  add("component_covariance",c("Income","Prediction","$V_r$","$V_u$","$V_v$","$2C_{ru}$","$2C_{rv}$","$2C_{uv}$","Total"),
      unlist(covrows),189L,columns="@{}llrrrrrrr@{}")
  upper <- read("upper_dominant_case_components")
  rows <- vapply(seq_len(nrow(upper)),function(i) {
    z <- upper[i,]; m <- paste0("dnn_",z$mechanism)
    line(c(labels[m],number(z$e_hat,3),number(z$m0_hat/1000),number(z$m1_hat/1000),
      number(z$contrast/1000),number(z$ineligible_residual/1000),number(z$phi/1000)))
  },character(1))
  add("upper_case",c("Prediction","Propensity","Control mean","Eligible mean","Contrast","Control correction","Signal"),
      rows,3L,columns="@{}lrrrrrr@{}")
  for(part in seq_len(3L)) {
    ids <- list(profiles[1:3],profiles[4:5],profiles[6:7])[[part]]
    rows <- list()
    for(id in ids)for(s in grid) {
      a <- lapply(methods,function(m)one(orders,id,m,s)); z <- a[[2L]]
      if(averaged_ties)stopifnot(z$nearest_weight>0,z$nearest_weight<=s/accepted$n+1e-12) else
        close(z$nearest_weight,s/accepted$n)
      rows[[length(rows)+1L]] <- line(c(plabel(id),number(s,0),number(z$ess,1),
        vapply(a,function(r)estimate(r$tau,r$se),character(1)),number(z$radius_90_z,4),number(100*z$jk_top1_sq_share,1)))
    }
    add(paste0("localization_",part),c("Income","$s$","$N$","All-tree","Single-fit/OOB","Five-fold","Radius","Largest (\\%)"),
        unlist(rows),length(rows)*3L,columns="@{}lrrrrrrr@{}")
  }
  rows <- vapply(profiles,function(id) {
    a <- one(asinh,id,"dnn_fullsample"); b <- one(asinh,id,"dnn_crossfit")
    line(c(plabel(id),paste0("\\(",value(a$tau),"\\;(",value(a$se),")\\)"),interval(a$ci_lo,a$ci_hi),
      paste0("\\(",value(b$tau),"\\;(",value(b$se),")\\)"),interval(b$ci_lo,b$ci_hi)))
  },character(1))
  add("asinh_historical",c("Income","Single-fit/OOB","95\\% interval","Five-fold","95\\% interval"),rows,14L)
  selected <- estimates[estimates$profile_id%in%display,]
  stopifnot(all(selected$ci_lo<0&selected$ci_hi>0))
  core_orders <- orders[orders$profile_id%in%display&orders$method=="dnn_oob"&
    orders$s%in%c(397L,99L),]
  stopifnot(nrow(selected)==9L,nrow(core_orders)==6L)
  rows <- list()
  for(id in display)for(m in methods) {
    e <- one(selected,id,m); z <- one(info,id,m)
    rows[[length(rows)+1L]] <- line(c(plabel(id),labels[m],number(e$tau/1000),
      number(e$se/1000),number(100*z$jk_top1_sq_share,1)))
  }
  add("core_estimates",c("Income","Prediction","Estimate","SE","Largest (\\%)"),
    unlist(rows),9L,columns="@{}llrrr@{}")
  rows <- list()
  for(id in display)for(s in c(397L,99L)) {
    z <- one(core_orders,id,"dnn_oob",s)
    if(averaged_ties)stopifnot(z$nearest_weight>0,z$nearest_weight<=s/accepted$n+1e-12) else
      close(z$nearest_weight,s/accepted$n)
    rows[[length(rows)+1L]] <- line(c(plabel(id),number(s,0),number(z$tau/1000),
      number(z$se/1000),number(z$ess,1),number(z$radius_90_z,4),
      number(100*z$jk_top1_sq_share,1)))
  }
  add("core_localization",c("Income","$s$","Estimate","SE","$N$","Radius","Largest (\\%)"),
    unlist(rows),6L,columns="@{}lrrrrrr@{}")
  anchors <- list(n=accepted$n,default_s=397L,N=unique(estimates$ess),
    all_nine_default_interior_intervals_include_zero=TRUE,
    upper_information=info[info$profile_id==profiles[5],],upper_components=upper,
    upper_orders=orders[orders$profile_id==profiles[5]&orders$method=="dnn_oob",],
    lower_orders=orders[orders$profile_id==profiles[3]&orders$method=="dnn_oob",],
    positive_exploratory_OOB_intervals=orders[orders$method=="dnn_oob"&orders$s<397L&
      orders$profile_id%in%display&orders$tau-qnorm(.975)*orders$se>0,c("profile_id","s","tau","se")])
  encoded <- function(x)paste0(jsonlite::toJSON(x,auto_unbox=TRUE,pretty=TRUE,digits=NA,na="null"),"\n")
  output_text <- c(setNames(lapply(tables,function(x)paste0(paste(x,collapse="\n"),"\n")),
    vapply(names(tables),function(x)file.path(table_dir,x),character(1))),
    setNames(list(encoded(anchors)),file.path(output_dir,"prose_anchors.json")))
  figure_names <- c("raw_dollars_intervals.pdf","raw_dollars_concentration.pdf","raw_dollars_localization.pdf")
  # Redraw the manuscript copies from portable summaries; historical PDFs stay intact.
  profile_titles <- c("Income $19,413","Income $31,476","Income $48,584")
  figure_records <- setNames(lapply(figure_names[1:2],function(name)list(
    input=relative(file.path(evidence,if(name==figure_names[1])"estimates_all.csv" else "information.csv")),
    reference_figure=relative(file.path(evidence,name)),
    output=relative(file.path(output_dir,name)),profiles=display,device="cairo_pdf",
    dimensions_inches=c(8,4.5),income_labels=sub("Income ","",profile_titles,fixed=TRUE),
    n_folds=fit$config$n_folds,
    records=if(name==figure_names[1])selected[,c("profile_id","method","tau","se","ci_lo","ci_hi")] else
      info[info$profile_id%in%display,c("profile_id","method","jk_top1_sq_share","jk_top5_sq_share")]
  )),figure_names[1:2])
  localization_name <- "raw_dollars_localization.pdf"
  localization_path <- file.path(output_dir,localization_name)
  localization_rows <- orders[orders$profile_id%in%display&orders$method=="dnn_oob",]
  figure_records[[localization_name]] <- list(
    input=relative(file.path(evidence,"localization_all.csv")),
    reference_figure=relative(file.path(evidence,localization_name)),
    output=relative(localization_path),profiles=display,method="dnn_oob",
    reference_order=397L,alpha=fit$config$alpha,profile_titles=profile_titles,
    device="cairo_pdf",dimensions_inches=c(10.5,4.4),
    records=localization_rows[,c("profile_id","method","s","tau","se")])
  if(check_only) {
    stored <- json(file.path(output_dir,"manifest.json"))
    figure_hashes <- setNames(lapply(figure_names,function(name) {
      h <- stored$output_hashes[[relative(file.path(output_dir,name))]]
      if(is.null(h))stop("Missing figure hash: ",name,call.=FALSE)
      h
    }),figure_names)
  } else {
    temporary_figures <- tempfile("sipp-manuscript-")
    dir.create(temporary_figures)
    on.exit(unlink(temporary_figures,recursive=TRUE),add=TRUE)
    plot_sipp_401k_affine_summaries(estimates,info,temporary_figures,display,
      sub("Income ","",profile_titles,fixed=TRUE),fit$config$n_folds,"raw_dollars")
    plot_sipp_401k_localization(localization_rows,display,profile_titles,397L,
      fit$config$alpha,file.path(temporary_figures,localization_name))
    figure_sources <- setNames(file.path(temporary_figures,figure_names),figure_names)
    figure_hashes <- lapply(figure_sources,sha)
  }
  output_hashes <- lapply(output_text,function(x)digest::digest(x,serialize=FALSE,algo="sha256"))
  names(output_hashes) <- vapply(names(output_text),relative,character(1))
  for(name in figure_names)output_hashes[[relative(file.path(output_dir,name))]] <- figure_hashes[[name]]
  sources <- c(.sipp_manuscript_source,.sipp_localization_source,.sipp_affine_source,
    file.path(repo_root,"Codes/Applications/SIPP_401k/07_render_affine_manuscript.R"))
  source_hashes <- setNames(lapply(sources,sha),vapply(sources,relative,character(1)))
  manifest <- list(schema_version=1L,as_of="2026-10-13",run_id=accepted$run_id,
    fit_fingerprint=accepted$fit_fingerprint,estimation_performed=FALSE,simulation_performed=FALSE,
    native_RDS_read=FALSE,renderer="Codes/Applications/SIPP_401k/07_render_affine_manuscript.R",
    renderer_sources=source_hashes,selection=list(all_profiles=profiles,display_profiles=display,
      methods=methods,default_s=397L,declared_orders=c(298L,397L,496L,745L),
      exploratory_orders=c(99L,149L,198L,248L),localization_observations=168L,
      core_default=list(profiles=display,methods=methods,order=397L,records=selected),
      core_localization=list(profiles=display,method="dnn_oob",orders=c(397L,99L),records=core_orders),
      historical_asinh_methods=c("single-fit/OOB","five-fold")),
    interpretation="Pointwise normal-reference intervals and frozen fitted-array attribution; no causal, bias, coverage or order-selection certification.",
    input_hashes=inputs,output_hashes=output_hashes,tables=table_records,figures=figure_records)
  output_text[[file.path(output_dir,"manifest.json")]] <- encoded(manifest)
  if(check_only) {
    for(path in names(output_text)) {
      if(!file.exists(path)||!identical(paste0(paste(readLines(path,warn=FALSE),collapse="\n"),"\n"),output_text[[path]]))
        stop("Generated output differs or is stale: ",relative(path),call.=FALSE)
    }
    for(name in figure_names) {
      path <- file.path(output_dir,name)
      if(!file.exists(path)||sha(path)!=output_hashes[[relative(path)]])
        stop("Generated figure hash differs or is stale: ",relative(path),call.=FALSE)
    }
  } else {
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
    for(path in names(output_text))writeLines(output_text[[path]],path,useBytes=TRUE,sep="")
    for(name in figure_names)if(!file.copy(figure_sources[[name]],file.path(output_dir,name),overwrite=TRUE))
      stop("Cannot publish figure: ",name,call.=FALSE)
  }
  invisible(manifest)
}
