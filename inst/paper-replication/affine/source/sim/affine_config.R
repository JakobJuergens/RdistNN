# Deterministic settings for the affine study; independent of historical runners.
affine_integer_root <- function(value, power) {
  out <- floor(exp(log(value) / power))
  while ((out + 1)^power <= value) out <- out + 1
  while (out^power > value) out <- out - 1
  as.integer(out)
}

affine_model_profile <- function(design_id) {
  if(length(design_id)!=1L || !design_id%in%c("C0","C1","C2","C0H","C2S")) stop("Unknown design_id")
  list(base=switch(design_id,C0H="C0",C2S="C2",design_id),
       sparse=design_id=="C2S",heteroskedastic=design_id=="C0H")
}

affine_schedule <- function(n, design_id = "C0", tuning = list()) {
  if (length(n) != 1L || !is.finite(n) || n < 2 || n != floor(n)) stop("n must be an integer >=2")
  profile <- affine_model_profile(design_id)
  allowed <- c("M", "q_dim", "s", "order_mult", "lambda", "lambda_mult")
  if (length(setdiff(names(tuning), allowed))) stop("Unknown tuning option")
  M <- if (is.null(tuning$M)) max(1L, affine_integer_root(n, 12)) else tuning$M
  q_dim <- if (is.null(tuning$q_dim)) if (profile$base == "C2") 5L else 3L else tuning$q_dim
  mult <- if (is.null(tuning$order_mult)) 1 else tuning$order_mult
  s <- if (!is.null(tuning$s)) tuning$s else affine_integer_root(mult^4 * n^3, 4)
  lm <- if (is.null(tuning$lambda_mult)) 1 else tuning$lambda_mult
  lambda <- if (is.null(tuning$lambda)) lm / sqrt(n) else tuning$lambda
  if (!is.finite(M) || M < 1 || M != floor(M) || !q_dim %in% c(3,5) ||
      !is.finite(mult) || mult <= 0 || !is.finite(s) || s < 1 || s >= n || s != floor(s) ||
      !is.finite(lambda) || lambda <= 0) stop("Invalid deterministic schedule (requires 1<=s<n, lambda>0)")
  list(n = as.integer(n), M = as.integer(M), J = as.integer(M^4),
       q_dim = as.integer(q_dim), K = as.integer(q_dim * M^4), s = as.integer(s),
       lambda = lambda, order_mult = mult, propensity_bounds = c(.25,.75),
       outcome_bounds = c(-1,1), x = c(.5,.5), norm_id = "euclidean_original",
       design_id = design_id, coverage = if (design_id == "C0" && q_dim == 3 &&
         M==max(1L,affine_integer_root(n,12)) && lambda==1/sqrt(n) && s==affine_integer_root(n^3,4))
         "primitive_corollary" else if (design_id == "C0" && q_dim == 3)
         "tuning_diagnostic_not_named_corollary" else if (profile$base == "C2" && q_dim == 3)
         "misspecified_diagnostic" else "explicit_extension_not_named_corollary")
}

affine_config <- function(stage = "preflight", overrides = list()) {
  stage <- match.arg(stage, c("preflight", "timing", "pilot", "production"))
  study <- if(is.null(overrides$study)) "dnn" else match.arg(overrides$study,c("dnn","kernel_comparison","settings_comparison"))
  make <- function(id, ns, orders = 1) {
    z <- expand.grid(n = ns, order_mult = orders, KEEP.OUT.ATTRS = FALSE)
    z$design_id <- id; z
  }
  cells <- if(study=="settings_comparison") {
    ns <- if(stage=="preflight") c(64,128) else c(2048,8192)
    rbind(make("C0H",ns),make("C2S",ns))
  } else if(study=="kernel_comparison") {
    ns <- if(stage=="preflight") c(64,128) else c(512,2048,8192)
    rbind(make("C0",ns),make("C2",ns,c(.4,1)))
  } else if (stage == "preflight") rbind(make("C0",c(64,128)), make("C1",c(64,128)),
    make("C2",c(64,128),c(.4,1))) else rbind(make("C0",c(512,2048,4096,16384,65536)),
    make("C1",c(2048,4096,16384)), make("C2",c(512,2048,8192),c(.4,1)))
  reps <- switch(stage, preflight = 20L, timing = 13L, pilot = 200L, production = 10000L)
  cfg <- list(stage = stage, study=study,
      localization_methods=if(study=="dnn") "dnn" else c("dnn","kernel_gaussian"),
      bandwidth_rule=if(study=="settings_comparison") "population_matched_dnn_ess_design_density" else "population_matched_dnn_ess_uniform_square",
      cells = cells, replications = reps, seed_capacity = 20000L,
      master_seed = 20261011L, rng_kind = c("Mersenne-Twister","Inversion","Rejection"),
      alpha = .05, folds = 5L, workers = 1L, chunk_size = 5L,
      control_reps = if(study!="dnn") 0L else if (stage == "preflight") 5L else 1000L,
      intensive_reps = if (stage == "preflight") 3L else 200L,
      cpu_budget = switch(stage,preflight=300,timing=7200,pilot=21600,production=NA_real_),
      bootstrap_replications = 199L,
      schema_version = "affine-v2", campaign_id = paste0("affine_",if(study=="kernel_comparison") "kernel_" else if(study=="settings_comparison") "settings_" else "",stage))
  if (length(setdiff(names(overrides),names(cfg)))) stop("Unknown config override")
  for (key in names(overrides)) cfg[[key]] <- overrides[[key]]
  expected_bandwidth <- if(study=="settings_comparison" || any(cfg$cells$design_id=="C2S"))
    "population_matched_dnn_ess_design_density" else "population_matched_dnn_ess_uniform_square"
  if(!"bandwidth_rule"%in%names(overrides)) cfg$bandwidth_rule <- expected_bandwidth
  for (key in c("replications","seed_capacity","folds","workers","chunk_size")) {
    v <- cfg[[key]]
    if (length(v)!=1 || !is.finite(v) || v<1 || v!=floor(v)) stop(paste("Invalid",key))
  }
  if (cfg$replications > cfg$seed_capacity || cfg$folds != 5L || cfg$alpha != .05)
    stop("Invalid replication/seed capacity or fixed inference contract")
  if (!identical(cfg$rng_kind,c("Mersenne-Twister","Inversion","Rejection")) ||
      cfg$schema_version!="affine-v2" || cfg$stage!=stage || cfg$study!=study ||
      !identical(cfg$localization_methods,if(study=="dnn") "dnn" else c("dnn","kernel_gaussian")) ||
      cfg$bandwidth_rule!=expected_bandwidth) stop("Unsupported RNG/schema/study override")
  for (key in c("control_reps","intensive_reps","bootstrap_replications")) {
    v <- cfg[[key]]
    if(length(v)!=1 || !is.finite(v) || v<0 || v!=floor(v)) stop(paste("Invalid",key))
  }
  affine_seed_map(cfg$seed_capacity,cfg$master_seed,stage)
  cells <- cfg$cells
  required <- c("n","order_mult","design_id")
  derived <- c("J","K","coverage","cell_id","N_reference","kernel_h")
  allowed_cells <- c(required,"M","q_dim","s","lambda","lambda_mult",derived,"replications","control_replications")
  if (!is.data.frame(cells) || !all(required %in% names(cells)) || nrow(cells)<1 ||
      length(setdiff(names(cells),allowed_cells))) stop("Invalid or unknown cell settings")
  tuning_names <- intersect(c("M","q_dim","s","lambda","lambda_mult","order_mult"),names(cells))
  schedules <- lapply(seq_len(nrow(cells)),function(i) {
    tuning <- as.list(cells[i,tuning_names,drop=FALSE])
    tuning <- tuning[!vapply(tuning,function(v) length(v)!=1 || is.na(v),logical(1))]
    affine_schedule(cells$n[i],cells$design_id[i],tuning)
  })
  for (key in c("M","J","q_dim","K","s","lambda","coverage")) {
    values <- vapply(schedules,`[[`,schedules[[1]][[key]],key)
    if(key%in%c("J","K","coverage") && key%in%names(cells) && !all(cells[[key]]==values))
      stop(paste("Inconsistent derived cell setting",key))
    cells[[key]] <- values
  }
  reference <- vapply(seq_len(nrow(cells)),function(i) 1/sum(affine_rank_weights(cells$n[i],cells$s[i])^2),numeric(1))
  h <- if(study!="dnn") vapply(seq_len(nrow(cells)),function(i)
      affine_kernel_bandwidth(cells$n[i],cells$s[i],design_id=cells$design_id[i]),numeric(1)) else rep(NA_real_,nrow(cells))
  for(key in c("N_reference","kernel_h")) {
    values <- if(key=="N_reference") reference else h
    if(key%in%names(cells) && !identical(as.numeric(cells[[key]]),values)) stop(paste("Inconsistent derived cell setting",key))
    cells[[key]] <- values
  }
  num_id <- function(v) format(v,digits=17,scientific=FALSE,trim=TRUE)
  ids <- paste(cells$design_id,cells$n,paste0("M",cells$M),paste0("q",cells$q_dim),
      paste0("s",cells$s),paste0("l",num_id(cells$lambda)),sep="_")
  if("cell_id"%in%names(cells) && !identical(cells$cell_id,ids)) stop("Cell IDs must match resolved tuning")
  cells$cell_id <- ids
  if (anyDuplicated(cells$cell_id)) stop("Duplicate cell identifiers")
  if(!"replications"%in%names(cells)) cells$replications <- cfg$replications
  if(!"control_replications"%in%names(cells)) cells$control_replications <- pmin(cfg$control_reps,cells$replications)
  if(any(!is.finite(cells$replications) | cells$replications<1 | cells$replications!=floor(cells$replications) |
      cells$replications>cfg$seed_capacity) || any(!is.finite(cells$control_replications) |
      cells$control_replications<0 | cells$control_replications!=floor(cells$control_replications) |
      cells$control_replications>cells$replications)) stop("Invalid per-cell replication budgets")
  cfg$cells <- cells; cfg

}

affine_preserve_rng <- function(expr) {
  present <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  if (present) old_seed <- get(".Random.seed", .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({do.call(RNGkind,as.list(old_kind)); if (present) assign(".Random.seed",old_seed,.GlobalEnv)
    else if (exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)}, add=TRUE)
  force(expr)
}

affine_with_seed <- function(seed,expr) affine_preserve_rng({
  RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(seed)
  force(expr)
})

affine_seed_map <- function(capacity = 20000L, master_seed = 20261011L,stage="production") {
  if (capacity < 1 || capacity != floor(capacity) || 3*capacity > 1e8) stop("Invalid seed capacity")
  stages <- c("preflight","timing","pilot","production")
  if(!stage%in%stages || 12*capacity>1e8) stop("Invalid reserved seed map")
  pool <- affine_with_seed(master_seed,matrix(sample.int(100000000L,12L*capacity,replace=FALSE),4L*capacity,3,byrow=TRUE))
  block <- (match(stage,stages)-1L)*capacity+seq_len(capacity)
  seeds <- pool[block,,drop=FALSE]
  data.frame(replication_id=seq_len(capacity),data_seed=seeds[,1],fold_seed=seeds[,2],integration_seed=seeds[,3])
}
