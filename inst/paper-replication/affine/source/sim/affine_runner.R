# Atomic run-local receipts. Numerical failures are rows; incomplete jobs are not.
affine_load_dependency <- function() {
  path <- Sys.getenv("SIM_RDISTNN_SOURCE_DIR",unset="")
  if(nzchar(path)) {
    path <- normalizePath(path,mustWork=TRUE)
    if(!identical(getOption("affine.dependency.source"),path)) {
      pkgload::load_all(path,quiet=TRUE,export_all=FALSE,helpers=FALSE)
      options(affine.dependency.source=path)
    }
  } else if(!requireNamespace("RdistNN",quietly=TRUE)) stop("RdistNN is required")
  if(!"jackknife_rank_full"%in%getNamespaceExports("RdistNN")) stop("RdistNN lacks its public full-rank covariance engine")
  invisible(TRUE)
}

affine_schedule_cell <- function(cell) {
  affine_schedule(cell$n,cell$design_id,list(M=cell$M,q_dim=cell$q_dim,s=cell$s,
      lambda=cell$lambda,order_mult=cell$order_mult))
}

affine_diagnostic_key <- function(cell_id,replication_id,method,localization_method="dnn")
  paste(cell_id,replication_id,method,localization_method,sep="/")

affine_base_row <- function(cfg,cell,replication_id,seeds,method,application,localization_method="dnn",local=NULL) {
  row <- data.frame(campaign_id=cfg$campaign_id,cell_id=cell$cell_id,replication_id=replication_id,
    application=application,training_method=method,localization_method=localization_method,
    data_seed=seeds$data_seed,fold_seed=seeds$fold_seed,
    integration_seed=seeds$integration_seed,design_id=cell$design_id,n=cell$n,M=cell$M,J=cell$J,q_dim=cell$q_dim,
    K=cell$K,s=cell$s,lambda=cell$lambda,x1=.5,x2=.5,
    coverage_status=if(localization_method=="dnn") cell$coverage else "kernel_extension_not_named_dnn_corollary",
    learner_coverage_status=cell$coverage,
    kernel_h=if(localization_method=="dnn") NA_real_ else cell$kernel_h,
    bandwidth_rule=if(localization_method=="dnn") "dnn_order" else cfg$bandwidth_rule,
    schema_version=cfg$schema_version,point_target=1/8,estimate=NA_real_,Gamma_hat=NA_real_,
    N=if(localization_method=="dnn") cell$N_reference else if(!is.null(local) && !inherits(local,"error")) local$N else NA_real_,
    N_reference=cell$N_reference,weighted_radius=NA_real_,largest_weight=NA_real_,score_mean_residual=NA_real_,
    weighted_truth=NA_real_,spatial_drift=NA_real_,centered_remainder=NA_real_,
    diagnostic_status="unavailable",diagnostic_error=NA_character_,
    status="failed",error_category=NA_character_,error_message=NA_character_,
    error_component=NA_character_,error_block=NA_integer_,error_fold=NA_integer_,gram_event=NA,
    min_raw_gram=NA_real_,empty_arm_cells=NA_integer_,max_clip_fraction=NA_real_)
  for(route in c("jackknife","direct_moment","direct_squared")) {
    for(field in c("variance","SE","lower","upper","score_cov_scaled")) row[[paste0(field,"_",route)]] <- NA_real_
    row[[paste0("variance_status_",route)]] <- "unavailable"
    row[[paste0("variance_error_",route)]] <- NA_character_
  }
  row
}

affine_method_rows <- function(cfg,cell,replication_id,seeds,method,compute,localization_method="dnn",local=NULL) {
  rows <- lapply(c("cate","slope"),function(app) affine_base_row(cfg,cell,replication_id,seeds,method,app,localization_method,local))
  state <- tryCatch(compute(),error=identity)
  if(inherits(state,"error")) {
    message <- conditionMessage(state)
    category <- if(grepl("propensity|prediction|nonfinite",message)) "invalid_predictions" else
      if(grepl("singular|positive definite|chol",message)) "solve_failure" else "method_error"
    for(i in seq_along(rows)) {
      rows[[i]]$error_category <- category; rows[[i]]$error_message <- message
      if(inherits(state,"affine_component_error")) {
        rows[[i]]$error_component <- state$component; rows[[i]]$error_block <- state$block; rows[[i]]$error_fold <- state$fold
      }
    }
    return(list(raw=do.call(rbind,rows),diagnostics=list(error=message),timing=data.frame()))
  }
  for(i in seq_along(rows)) {
    app <- c("cate","slope")[i]; a <- state$inference[[app]]; r <- rows[[i]]
    r$N <- state$local$N; r$weighted_radius <- state$local$weighted_radius; r$largest_weight <- state$local$largest_weight
    r$estimate <- a$estimate; r$Gamma_hat <- a$Gamma_hat; r$status <- a$point_status
    r$score_mean_residual <- sum(state$local$weights*a$score)
    d <- state$diagnostics
    r$diagnostic_status <- if(!is.null(d$error)) "failed" else "ok"
    r$diagnostic_error <- if(!is.null(d$error)) d$error else NA_character_
    r$weighted_truth <- state$weighted_truth[app]; r$spatial_drift <- r$weighted_truth-r$point_target
    r$centered_remainder <- a$estimate-r$weighted_truth
    if(!is.null(state$fit$diagnostics)) {
      r$gram_event <- state$fit$diagnostics$gram_event
      r$min_raw_gram <- min(state$fit$diagnostics$min_raw)
      r$empty_arm_cells <- sum(state$fit$diagnostics$empty_counts[c("mu0","mu1")])
    }
    if(!is.null(state$predictions$clip_fraction)) r$max_clip_fraction <- max(state$predictions$clip_fraction)
    for(route in names(a$variance)) {
      r[[paste0("variance_",route)]] <- unname(a$variance[route]); r[[paste0("SE_",route)]] <- unname(a$SE[route])
      r[[paste0("lower_",route)]] <- a$intervals[route,1]; r[[paste0("upper_",route)]] <- a$intervals[route,2]
      r[[paste0("score_cov_scaled_",route)]] <- unname(a$score_cov_scaled[route])
      r[[paste0("variance_status_",route)]] <- unname(a$variance_status[route])
      r[[paste0("variance_error_",route)]] <- unname(a$variance_error[route])
    }
    rows[[i]] <- r
  }
  list(raw=do.call(rbind,rows),diagnostics=state$diagnostics,timing=state$timing)
}

affine_family <- function(replication_id,cfg,seed_map) {
  # Rcpp's RNGScope may initialize a seed even for a deterministic engine call.
  affine_preserve_rng(affine_family_impl(replication_id,cfg,seed_map))
}

affine_family_impl <- function(replication_id,cfg,seed_map) {
  affine_load_dependency(); start_cpu <- sum(proc.time()[1:2]); start_wall <- proc.time()[3]
  seeds <- seed_map[replication_id,,drop=FALSE]
  latent <- affine_latent(max(cfg$cells$n),seeds$data_seed)
  weights_cache <- fit_cache <- new.env(parent=emptyenv())
  # Separate environments prevent accidental dictionary/localizer key collisions.
  weights_cache <- new.env(parent=emptyenv())
  raw <- diag <- times <- list()
  for(ci in seq_len(nrow(cfg$cells))) {
    cell <- cfg$cells[ci,,drop=FALSE]; sch <- affine_schedule_cell(cell)
    if(replication_id>cell$replications) next
    data <- affine_data(latent[seq_len(cell$n),,drop=FALSE],cell$design_id)
    for(localization_method in cfg$localization_methods) {
    geometry <- if(affine_model_profile(cell$design_id)$sparse) "sparse" else "uniform"
    local_key <- paste(geometry,localization_method,cell$n,cell$s,cell$kernel_h,sep="_")
    local_start <- proc.time()[3]; local_reused <- exists(local_key,weights_cache,inherits=FALSE)
    if(!local_reused) assign(local_key,tryCatch(if(localization_method=="dnn") affine_weights(data$X,sch$x,sch$s)
        else affine_kernel_weights(data$X,sch$x,cell$kernel_h,sch$s,cell$design_id),error=identity),weights_cache)
    local <- get(local_key,weights_cache); local_time <- unname(proc.time()[3]-local_start)
    methods <- c("full","pooled_five_fold","oracle")
    if(replication_id<=cell$control_replications) methods <- c(methods,"population_projection","population_ridge")
    for(method in methods) {
      compute <- function() {
        if(inherits(local,"error")) stop(local)
        key <- paste(cell$design_id,cell$n,cell$M,cell$q_dim,cell$lambda,method,sep="_")
        reused <- exists(key,fit_cache,inherits=FALSE); begin <- proc.time()[3]
        if(!reused) {
          des <- affine_design(data$X,data$Z,sch)
          value <- switch(method,
            full={f <- affine_fit(des,data$T,data$Y,timing=TRUE); list(fit=f,predictions=affine_predict(f,data$X,data$Z))},
            pooled_five_fold={full_key <- sub("pooled_five_fold$","full",key)
              full_stats <- if(exists(full_key,fit_cache,inherits=FALSE)) get(full_key,fit_cache)$fit$stats else NULL
              f <- affine_crossfit(des,data$T,data$Y,seeds$fold_seed,full_stats=full_stats,timing=TRUE)
              list(fit=f,predictions=f$predictions)},
            oracle=list(fit=list(method="oracle",schedule=sch),predictions=data$truth),
            population_projection={f <- affine_population_fit(sch,cell$design_id); list(fit=f,predictions=affine_predict(f,data$X,data$Z))},
            population_ridge={f <- affine_population_fit(sch,cell$design_id,ridge=TRUE); list(fit=f,predictions=affine_predict(f,data$X,data$Z))})
          assign(key,value,fit_cache)
        }
        value <- get(key,fit_cache); first_time <- unname(proc.time()[3]-begin)
        begin <- proc.time()[3]; eq <- affine_equations(data,value$predictions); equation_time <- unname(proc.time()[3]-begin)
        begin <- proc.time()[3]; inf <- affine_inference(eq,local); covariance_time <- unname(proc.time()[3]-begin)
        intensive_cell <- cfg$stage=="preflight" || cell$design_id%in%c("C0H","C2S") || (cell$design_id=="C0" && cell$n%in%c(2048,4096,65536)) ||
             (cell$design_id=="C1" && cell$n==4096) || (cell$design_id=="C2" && cell$n==2048)
        intensive <- replication_id<=cfg$intensive_reps && intensive_cell && method%in%c("full","pooled_five_fold","oracle")
        begin <- proc.time()[3]
        v0 <- data$T-data$truth$pi
        weighted_truth <- c(cate=sum(local$weights*data$truth$theta),
          slope=sum(local$weights*v0^2*data$truth$theta)/sum(local$weights*v0^2))
        d <- tryCatch(affine_diagnostics(data,value$fit,eq,local,cell$design_id,
          list(intensive=intensive),inference=inf),error=function(e) list(error=conditionMessage(e)))
        diagnostic_time <- unname(proc.time()[3]-begin)
        timing <- data.frame(replication_id=replication_id,cell_id=cell$cell_id,training_method=method,localization_method=localization_method,
          first_stage=first_time,equations=equation_time,localization=local_time,covariance=covariance_time,
          diagnostics=diagnostic_time,weights_reused=local_reused,first_stage_reused=reused,backend=if(cfg$workers==1) "serial" else "multisession")
        component_seconds <- if(reused) setNames(numeric(4),c("pi","mu0","mu1","m")) else
          if(!is.null(value$fit$fits)) Reduce(`+`,lapply(value$fit$fits,`[[`,"component_fit_seconds")) else value$fit$component_fit_seconds
        if(is.null(component_seconds)) component_seconds <- setNames(numeric(4),c("pi","mu0","mu1","m"))
        for(j in names(component_seconds)) timing[[paste0("solve_",j)]] <- component_seconds[j]
        list(fit=value$fit,predictions=value$predictions,inference=inf,local=local,weighted_truth=weighted_truth,diagnostics=d,timing=timing)
      }
      result <- affine_method_rows(cfg,cell,replication_id,seeds,method,compute,localization_method,local)
      raw[[length(raw)+1L]] <- result$raw
      diag[[affine_diagnostic_key(cell$cell_id,replication_id,method,localization_method)]] <- result$diagnostics
      if(nrow(result$timing)) times[[length(times)+1L]] <- result$timing
    }
    }
  }
  list(replication_id=replication_id,raw=do.call(rbind,raw),diagnostics=diag,
       timing=if(length(times)) do.call(rbind,times) else data.frame(),
       worker_cpu_seconds=sum(proc.time()[1:2])-start_cpu,worker_wall_seconds=unname(proc.time()[3]-start_wall))
}

affine_expected_keys <- function(cfg,ids) {
  out <- list()
  for(r in ids) for(i in seq_len(nrow(cfg$cells))) {
    c <- cfg$cells[i,]; if(r>c$replications) next
    methods <- c("full","pooled_five_fold","oracle")
    if(r<=c$control_replications) methods <- c(methods,"population_projection","population_ridge")
    z <- expand.grid(application=c("cate","slope"),training_method=methods,
                     localization_method=cfg$localization_methods,stringsAsFactors=FALSE)
    z$cell_id <- c$cell_id; z$replication_id <- r; out[[length(out)+1L]] <- z
  }
  do.call(rbind,out)
}

affine_row_keys <- function(raw) do.call(paste,c(raw[c("cell_id","replication_id","application","training_method","localization_method")],sep="/"))

affine_payload_hash <- function(chunk) digest::digest(chunk[setdiff(names(chunk),"payload_hash")],algo="sha256")

affine_raw_schema <- function(cfg,seed_map) {
  template <- affine_base_row(cfg,cfg$cells[1,,drop=FALSE],1L,seed_map[1,,drop=FALSE],"full","cate")
  type <- function(v) if(is.numeric(v)) "numeric" else typeof(v)
  fields <- lapply(template,function(v) list(type=type(v),nullable=TRUE))
  for(k in c("campaign_id","cell_id","replication_id","application","training_method","data_seed",
      "fold_seed","integration_seed","design_id","n","M","J","q_dim","K","s","lambda",
      "localization_method","bandwidth_rule","learner_coverage_status","N_reference",
      "x1","x2","coverage_status","schema_version","point_target","status","diagnostic_status"))
    fields[[k]]$nullable <- FALSE
  fields$status$values <- c("ok","failed","unavailable")
  fields$diagnostic_status$values <- c("ok","failed","unavailable")
  for(route in c("jackknife","direct_moment","direct_squared"))
    fields[[paste0("variance_status_",route)]]$values <- c("ok","zero","unavailable")
  fields$localization_method$values <- cfg$localization_methods
  list(version=cfg$schema_version,primary_key=c("campaign_id","cell_id","replication_id","application","training_method","localization_method"),
    raw_fields=fields,variance_methods=c("jackknife","direct_moment","direct_squared"),quantile_type=8L,
    definitions=list(N="realized 1/sum(w^2), nullable only on localization failure",
      N_reference="deterministic exact-DNN ESS; Gaussian population ESS is matched to this",
      kernel_h="deterministic common isotropic bandwidth; NA for DNN",
      direct_squared="analytic kernel sandwich for kernel_gaussian; diagnostic for dnn",
      point_target="theta(0.5,0.5)=1/8",
      variance="actual target variance",score_cov_scaled="N * Gamma_hat^2 * target variance",
      frozen_jackknife="fixed nuisances and target; same s for DNN, same h and renormalized retained weights for kernel",
      interval_denominator="successful estimate and finite available SE >= 0",
      student_denominator="available interval and SE > 0",missing="NA, accompanied by status/error",
      timing="separate payload; excluded from raw/summary numerical hashes",
      chunk="freeze hash, replication IDs, families with raw/diagnostics/timing, payload SHA256"))
}

affine_validate_raw <- function(raw,cfg,seed_map) {
  schema <- affine_raw_schema(cfg,seed_map)
  if(!is.data.frame(raw) || !identical(names(raw),names(schema$raw_fields))) stop("Raw schema fields mismatch")
  for(k in names(schema$raw_fields)) {
    f <- schema$raw_fields[[k]]; v <- raw[[k]]
    type_ok <- if(f$type=="numeric") is.numeric(v) else typeof(v)==f$type
    if(!type_ok || (!f$nullable && anyNA(v)) || (!is.null(f$values) && !all(v%in%f$values)))
      stop(paste("Invalid raw field",k))
  }
  cells <- cfg$cells[match(raw$cell_id,cfg$cells$cell_id),,drop=FALSE]
  seeds <- seed_map[match(raw$replication_id,seed_map$replication_id),,drop=FALSE]
  if(anyNA(cells$cell_id) || anyNA(seeds$replication_id)) stop("Unknown raw cell/ID")
  for(k in c("n","M","J","q_dim","K","s","lambda","design_id"))
    if(!all(raw[[k]]==cells[[k]])) stop(paste("Frozen tuning mismatch",k))
  for(k in c("data_seed","fold_seed","integration_seed"))
    if(!all(raw[[k]]==seeds[[k]])) stop(paste("Frozen seed mismatch",k))
  kernel <- raw$localization_method=="kernel_gaussian"
  expected_coverage <- ifelse(kernel,"kernel_extension_not_named_dnn_corollary",cells$coverage)
  if(!all(raw$campaign_id==cfg$campaign_id & raw$schema_version==cfg$schema_version &
      raw$coverage_status==expected_coverage & raw$learner_coverage_status==cells$coverage &
      raw$x1==.5 & raw$x2==.5 & raw$point_target==1/8 & raw$N_reference==cells$N_reference) ||
      any(!is.na(raw$kernel_h[!kernel])) || !all(raw$kernel_h[kernel]==cells$kernel_h[kernel]) ||
      !all(raw$bandwidth_rule==ifelse(kernel,cfg$bandwidth_rule,"dnn_order")))
    stop("Frozen reporting metadata mismatch")
  if(!identical(raw$N[!kernel],cells$N_reference[!kernel]) ||
     any(!is.finite(raw$N) & raw$status=="ok") ||
     any(raw$N<1-1e-10 | raw$N>raw$n+1e-10,na.rm=TRUE) ||
     any(raw$status=="ok" & !is.finite(raw$estimate))) stop("Invalid effective size/point status")
  # Reconstruct each distinct realized kernel ESS once, rather than accepting
  # a duplicated reference ESS or trusting a resealed numerical payload.
  groups <- split(which(kernel & !is.na(raw$N)),paste(raw$cell_id[kernel & !is.na(raw$N)],
     raw$replication_id[kernel & !is.na(raw$N)],sep="/"))
  for(ix in groups) {
    i <- ix[1]; X <- affine_reporting_quantile(affine_latent(raw$n[i],raw$data_seed[i])[,1:2,drop=FALSE],raw$design_id[i])
    expected_N <- affine_kernel_weights(X,h=raw$kernel_h[i],s=raw$s[i],design_id=raw$design_id[i])$N
    if(!all(raw$N[ix]==expected_N)) stop("Invalid realized kernel effective size")
  }
  invisible(TRUE)
}

affine_validate_chunk <- function(chunk,cfg,hash,seed_map=NULL) {
  if(!identical(chunk$freeze_hash,hash)) stop("Chunk freeze hash mismatch")
  if(!identical(chunk$payload_hash,affine_payload_hash(chunk))) stop("Chunk payload hash mismatch")
  ids <- chunk$replication_ids
  if(!is.integer(ids) || !length(ids) || anyNA(ids) || any(ids<1 | ids>max(cfg$cells$replications)) ||
      anyDuplicated(ids)) stop("Invalid chunk replication IDs")
  if(is.null(seed_map)) seed_map <- affine_seed_map(cfg$seed_capacity,cfg$master_seed,cfg$stage)
  affine_validate_raw(chunk$raw,cfg,seed_map)
  expected <- affine_expected_keys(cfg,ids); keys <- affine_row_keys(chunk$raw)
  if(anyDuplicated(keys) || !setequal(keys,affine_row_keys(expected))) stop("Incomplete or inconsistent chunk keys")
  if(!is.list(chunk$families) || length(chunk$families)!=length(ids) ||
      !identical(vapply(chunk$families,`[[`,integer(1),"replication_id"),ids)) stop("Chunk/family IDs mismatch")
  combined <- do.call(rbind,lapply(chunk$families,`[[`,"raw"))
  if(!identical(combined,chunk$raw)) stop("Family/raw payload mismatch")
  for(f in chunk$families) {
    if(!all(f$raw$replication_id==f$replication_id)) stop("Family raw IDs mismatch")
    expected_diagnostics <- unique(affine_diagnostic_key(f$raw$cell_id,f$replication_id,
       f$raw$training_method,f$raw$localization_method))
    if(!setequal(names(f$diagnostics),expected_diagnostics) || anyDuplicated(names(f$diagnostics)) ||
       !is.data.frame(f$timing) || !is.finite(f$worker_cpu_seconds) || f$worker_cpu_seconds<0 ||
       !is.finite(f$worker_wall_seconds) || f$worker_wall_seconds<0) stop("Missing family diagnostic/timing receipts")
  }
  invisible(TRUE)
}

affine_atomic_rds <- function(object,path) {
  temp <- paste0(path,".tmp-",Sys.getpid()); saveRDS(object,temp,version=3)
  if(!file.rename(temp,path)) stop("Atomic chunk rename failed")
}

affine_freeze <- function(cfg,source_dir,seed_map) {
  scientific <- cfg; scientific[c("workers","chunk_size","cpu_budget")] <- NULL
  files <- sort(list.files(source_dir,pattern="\\.R$",full.names=TRUE))
  source_hashes <- vapply(files,function(f) digest::digest(file=f,algo="sha256"),character(1))
  names(source_hashes) <- basename(files)
  package_dir <- find.package("RdistNN")
  source_override <- Sys.getenv("SIM_RDISTNN_SOURCE_DIR",unset="")
  if(nzchar(source_override)) {
    # Hash executable sources and metadata, not unrelated notes, figures or .git.
    package_files <- c(file.path(package_dir,c("DESCRIPTION","NAMESPACE")),
      list.files(file.path(package_dir,"R"),pattern="\\.R$",full.names=TRUE),
      list.files(file.path(package_dir,"src"),pattern="\\.(cpp|cc|c|h|hpp)$|^Makevars",full.names=TRUE))
  } else package_files <- list.files(package_dir,recursive=TRUE,full.names=TRUE)
  package_files <- sort(package_files[file.exists(package_files) & !dir.exists(package_files)])
  dll_paths <- vapply(getLoadedDLLs(),function(d) d[["path"]],character(1))
  dll_paths <- dll_paths[grepl("RdistNN",basename(dll_paths),fixed=TRUE)]
  file_hashes <- vapply(package_files,function(f) digest::digest(file=f,algo="sha256"),character(1))
  names(file_hashes) <- substring(package_files,nchar(package_dir)+2L)
  git_commit <- NA_character_
  if(nzchar(source_override) && nzchar(Sys.which("git"))) {
    git_commit <- suppressWarnings(system2("git",c("-C",shQuote(package_dir),"rev-parse","HEAD"),stdout=TRUE,stderr=FALSE))
    if(!is.null(attr(git_commit,"status"))) git_commit <- NA_character_
  }
  dependency <- list(version=as.character(packageVersion("RdistNN")),path=package_dir,files=file_hashes,
      loaded_DLL_hashes=vapply(dll_paths,function(f) digest::digest(file=f,algo="sha256"),character(1)),
      public_engine=digest::digest(body(RdistNN::jackknife_rank_full),algo="sha256"),
      source_override=source_override,source_commit=git_commit)
  freeze <- list(config=scientific,seed_hash=digest::digest(seed_map,algo="sha256"),
       source_hashes=source_hashes,dependency=dependency,R_version=R.version.string)
  list(hash=digest::digest(freeze,algo="sha256"),details=freeze)
}

affine_run <- function(cfg,output_dir,resume=FALSE,max_chunks=Inf,source_dir=NULL) {
  affine_load_dependency()
  invisible(affine_config(cfg$stage,cfg))
  if(is.null(source_dir)) source_dir <- file.path(sim_find_project_root(),"Codes","Simulations","affine_scores")
  for(i in seq_len(nrow(cfg$cells))) affine_schedule_cell(cfg$cells[i,])
  if(!is.finite(cfg$cpu_budget) || cfg$cpu_budget<=0) stop("Supply a positive accepted worker CPU budget before production")
  seed_map <- affine_seed_map(cfg$seed_capacity,cfg$master_seed,cfg$stage)
  freeze <- affine_freeze(cfg,source_dir,seed_map)
  if(dir.exists(output_dir) && !resume) stop("Output directory exists; use a new run name or explicit identical resume")
  if(resume) {
    old <- readRDS(file.path(output_dir,"freeze.rds"))
    if(!identical(old$hash,freeze$hash)) stop("Resume design/learner/seed/schema/source/dependency mismatch")
  } else {
    dir.create(file.path(output_dir,"chunks"),recursive=TRUE)
    saveRDS(freeze,file.path(output_dir,"freeze.rds")); saveRDS(cfg,file.path(output_dir,"config.rds"))
    jsonlite::write_json(cfg,file.path(output_dir,"config.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
    write.csv(cfg$cells,file.path(output_dir,"cells.csv"),row.names=FALSE)
    write.csv(seed_map,file.path(output_dir,"seed_map.csv"),row.names=FALSE)
    jsonlite::write_json(affine_raw_schema(cfg,seed_map),file.path(output_dir,"schema.json"),auto_unbox=TRUE,pretty=TRUE)
    writeLines(capture.output(sessionInfo()),file.path(output_dir,"session.txt"))
    writeLines(c("# Affine run receipt",paste("Stage:",cfg$stage),paste("Freeze hash:",freeze$hash),
      "Pilot/preflight rows are diagnostics and cannot be promoted as accepted production."),file.path(output_dir,"README.md"))
  }
  chunks <- lapply(list.files(file.path(output_dir,"chunks"),pattern="\\.rds$",full.names=TRUE),readRDS)
  for(ch in chunks) affine_validate_chunk(ch,cfg,freeze$hash,seed_map)
  completed <- unlist(lapply(chunks,`[[`,"replication_ids"))
  if(anyDuplicated(completed)) stop("Overlapping completed chunks")
  remaining <- setdiff(seq_len(max(cfg$cells$replications)),completed)
  jobs <- split(remaining,ceiling(seq_along(remaining)/cfg$chunk_size))
  total_cpu <- sum(vapply(chunks,function(ch) sum(vapply(ch$families,`[[`,numeric(1),"worker_cpu_seconds")),numeric(1)))
  old_plan <- future::plan(); on.exit(future::plan(old_plan),add=TRUE)
  old_threads <- Sys.getenv(c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"),unset=NA_character_)
  on.exit({for(key in names(old_threads)) if(is.na(old_threads[[key]])) Sys.unsetenv(key) else do.call(Sys.setenv,setNames(list(old_threads[[key]]),key))},add=TRUE)
  Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
  if(cfg$workers>1) future::plan(future::multisession,workers=cfg$workers) else future::plan(future::sequential)
  jsonlite::write_json(list(status="running",freeze_hash=freeze$hash,completed_replications=length(completed)),file.path(output_dir,"manifest.json"),auto_unbox=TRUE,pretty=TRUE)
  for(j in seq_along(jobs)) {
    if(j>max_chunks || total_cpu>=cfg$cpu_budget) break
    ids <- jobs[[j]]
    families <- if(cfg$workers==1) lapply(ids,affine_family,cfg=cfg,seed_map=seed_map) else
       future.apply::future_lapply(ids,affine_family,cfg=cfg,seed_map=seed_map,future.seed=FALSE)
    chunk <- list(freeze_hash=freeze$hash,replication_ids=as.integer(ids),families=families,
                  raw=do.call(rbind,lapply(families,`[[`,"raw")))
    chunk$payload_hash <- affine_payload_hash(chunk)
    affine_validate_chunk(chunk,cfg,freeze$hash,seed_map)
    affine_atomic_rds(chunk,file.path(output_dir,"chunks",sprintf("chunk-%06d-%06d.rds",min(ids),max(ids))))
    chunks[[length(chunks)+1L]] <- chunk
    total_cpu <- total_cpu+sum(vapply(families,`[[`,numeric(1),"worker_cpu_seconds"))
    completed_now <- sum(vapply(chunks,function(ch) length(ch$replication_ids),integer(1)))
    throughput <- if(total_cpu>0) completed_now/total_cpu else NA_real_
    jsonlite::write_json(list(status="running",freeze_hash=freeze$hash,completed_replications=completed_now,
      planned_replications=max(cfg$cells$replications),worker_cpu_seconds=total_cpu,
      families_per_worker_cpu_second=throughput,
      projected_remaining_worker_cpu_seconds=(max(cfg$cells$replications)-completed_now)/throughput),
      file.path(output_dir,"manifest.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
  }
  families <- unlist(lapply(chunks,`[[`,"families"),recursive=FALSE)
  raw <- if(length(chunks)) do.call(rbind,lapply(chunks,`[[`,"raw")) else data.frame()
  if(nrow(raw)) {raw <- raw[do.call(order,raw[c("cell_id","replication_id","application","training_method","localization_method")]),]; rownames(raw) <- NULL}
  diagnostics <- unlist(lapply(families,`[[`,"diagnostics"),recursive=FALSE)
  timing <- if(length(families)) do.call(rbind,lapply(families,`[[`,"timing")) else data.frame()
  completed <- sort(unlist(lapply(chunks,`[[`,"replication_ids")))
  status <- if(identical(as.integer(completed),seq_len(max(cfg$cells$replications)))) "completed" else "incomplete"
  summary <- affine_summarise(raw,diagnostics,cfg$cells,bootstrap=cfg$bootstrap_replications)
  saveRDS(raw,file.path(output_dir,"raw.rds")); saveRDS(diagnostics,file.path(output_dir,"diagnostics.rds")); saveRDS(timing,file.path(output_dir,"timing.rds"))
  saveRDS(summary,file.path(output_dir,"summary.rds"))
  write.csv(raw,file.path(output_dir,"raw.csv"),row.names=FALSE)
  for(key in names(summary)[vapply(summary,is.data.frame,logical(1))]) write.csv(summary[[key]],file.path(output_dir,paste0("summary_",key,".csv")),row.names=FALSE)
  manifest <- list(status=status,freeze_hash=freeze$hash,completed_replications=length(completed),
      planned_main_rows=sum(cfg$cells$replications)*6*length(cfg$localization_methods),
      planned_control_rows=sum(cfg$cells$control_replications)*4*length(cfg$localization_methods),
      actual_rows=nrow(raw),numerical_failures=if(nrow(raw)) sum(raw$status!="ok") else 0L,
      worker_cpu_seconds=total_cpu,cpu_budget=cfg$cpu_budget,budget_check="between_atomic_chunks",
      numerical_hash=digest::digest(raw,algo="sha256"),summary_hash=digest::digest(summary,algo="sha256"),
      diagnostic_failures=if(nrow(raw)) sum(raw$diagnostic_status=="failed") else 0L)
  schema <- affine_raw_schema(cfg,seed_map)
  schema$summary_fields <- lapply(summary[vapply(summary,is.data.frame,logical(1))],function(z)
      lapply(z,function(v) list(type=if(is.numeric(v)) "numeric" else typeof(v),nullable=TRUE)))
  schema$validation <- list(raw_fields_valid=if(nrow(raw)) isTRUE(affine_validate_raw(raw,cfg,seed_map)) else TRUE,
       raw_rows=nrow(raw),completed_chunk_count=length(chunks),validated_payload_hashes=TRUE)
  jsonlite::write_json(schema,file.path(output_dir,"schema.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
  jsonlite::write_json(manifest,file.path(output_dir,"manifest.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
  list(status=status,raw=raw,diagnostics=diagnostics,timing=timing,summary=summary,manifest=manifest)
}
