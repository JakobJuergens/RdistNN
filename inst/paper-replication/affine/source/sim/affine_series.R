# Sparse cell blocks implement the normalized dictionary objective exactly.
affine_design <- function(X,Z,schedule) {
  if (nrow(X)!=schedule$n || nrow(Z)!=schedule$n) stop("Design size differs from ambient n")
  list(Q=affine_basis(X,schedule$q_dim),cell=affine_cells(Z,schedule$M),schedule=schedule,n=schedule$n)
}

affine_validate_retained <- function(retain,n) {
  if (!is.numeric(retain) || anyNA(retain) || any(retain!=floor(retain)) ||
      any(retain<1 | retain>n) || anyDuplicated(retain)) stop("Invalid retained indices")
  as.integer(retain)
}

affine_stats <- function(design,T,Y,retain=seq_len(design$n)) {
  n <- design$n; qd <- design$schedule$q_dim; J <- design$schedule$J
  if (length(T)!=n || length(Y)!=n || anyNA(T) || any(!T%in%0:1) || any(!is.finite(Y))) stop("Invalid regression data")
  retain <- affine_validate_retained(retain,n)
  empty <- list(S=array(0,c(qd,qd,J)),h=matrix(0,qd,J),count=integer(J))
  out <- setNames(rep(list(empty),4),c("pi","mu0","mu1","m"))
  groups <- split(retain,design$cell[retain])
  for (cell_name in names(groups)) {
    c <- as.integer(cell_name); ix <- groups[[cell_name]]
    Q <- design$Q[ix,,drop=FALSE]; S <- crossprod(Q)
    out$pi$S[,,c] <- out$m$S[,,c] <- S
    out$pi$count[c] <- out$m$count[c] <- length(ix)
    out$pi$h[,c] <- crossprod(Q,T[ix]); out$m$h[,c] <- crossprod(Q,Y[ix])
    for (arm in 0:1) {
      a <- ix[T[ix]==arm]; j <- paste0("mu",arm)
      Qa <- design$Q[a,,drop=FALSE]
      out[[j]]$S[,,c] <- crossprod(Qa)
      out[[j]]$h[,c] <- crossprod(Qa,Y[a]); out[[j]]$count[c] <- length(a)
    }
  }
  out
}

affine_stats_difference <- function(full,removed) {
  lapply(names(full),function(j) list(S=full[[j]]$S-removed[[j]]$S,
     h=full[[j]]$h-removed[[j]]$h,count=full[[j]]$count-removed[[j]]$count)) |>
    setNames(names(full))
}

affine_solve_block <- function(A,h) {
  R <- chol(A)
  drop(backsolve(R,forwardsolve(t(R),h)))
}

affine_gram_diagnostics <- function(stats,schedule) {
  eigenvalues <- lapply(stats,function(st) vapply(seq_len(schedule$J),function(c)
    eigen(st$S[,,c]*(schedule$J/schedule$n),symmetric=TRUE,only.values=TRUE)$values,numeric(schedule$q_dim)))
  spectra <- lapply(eigenvalues,function(v) rbind(apply(v,2,min),apply(v,2,max)))
  min_raw <- vapply(spectra,function(v) min(v[1,]),numeric(1))
  # Tiny negative eigenvalues in singular empirical blocks are roundoff diagnostics.
  list(raw_eigenvalues=eigenvalues,penalized_eigenvalues=lapply(eigenvalues,function(v) v+schedule$lambda),
       min_raw=min_raw,min_penalized=min_raw+schedule$lambda,
       gram_event=all(min_raw>=3/16),empty_counts=vapply(stats,function(st) sum(st$count==0),integer(1)),
       max_raw_condition=vapply(spectra,function(v) if(any(v[1,]<=0)) Inf else max(v[2,]/v[1,]),numeric(1)),
       max_penalized_condition=vapply(spectra,function(v) max((v[2,]+schedule$lambda)/(v[1,]+schedule$lambda)),numeric(1)),
       two_deletion_perturbation=2*(if(schedule$q_dim==3) 7 else 11)*schedule$J/schedule$n)
}

affine_fit_from_stats <- function(stats,schedule,retain,timing=FALSE) {
  penalty <- diag(schedule$n*schedule$lambda/schedule$J,schedule$q_dim)
  component_seconds <- numeric(length(stats)) |> setNames(names(stats))
  coefs <- lapply(names(stats),function(j) {
    st <- stats[[j]]; begin <- if(timing) proc.time()[3] else 0
    c <- t(vapply(seq_len(schedule$J),function(c) tryCatch(
      affine_solve_block(st$S[,,c]+penalty,st$h[,c]),error=function(e) stop(structure(
        list(message=paste0("Regression ",j,", block ",c,": ",conditionMessage(e)),
          call=NULL,component=j,block=as.integer(c),fold=NA_integer_),
        class=c("affine_component_error","error","condition")))),numeric(schedule$q_dim)))
    if(timing) component_seconds[j] <<- unname(proc.time()[3]-begin)
    c
  }) |> setNames(names(stats))
  fit <- list(coefficients=coefs,schedule=schedule,stats=stats,retain=retain,
       diagnostics=affine_gram_diagnostics(stats,schedule),update_fallbacks=0L,
       coefficient_norms=vapply(coefs,function(c) sqrt(sum(c^2)/schedule$J),numeric(1)),
       parameterization="cell_q",method="sample_ridge")
  if(timing) fit$component_fit_seconds <- component_seconds
  fit
}

affine_fit <- function(design,T,Y,retain=seq_len(design$n),ambient_n=design$n,lambda=design$schedule$lambda,timing=FALSE) {
  if (ambient_n!=design$n || lambda!=design$schedule$lambda) stop("Retained fits must keep ambient n and lambda fixed")
  retain <- affine_validate_retained(retain,design$n)
  affine_fit_from_stats(affine_stats(design,T,Y,retain),design$schedule,retain,timing)
}

affine_predict <- function(fit,X,Z) {
  if (nrow(X)!=nrow(Z)) stop("Prediction row counts differ")
  Q <- affine_basis(X,fit$schedule$q_dim); cells <- affine_cells(Z,fit$schedule$M)
  raw <- lapply(fit$coefficients,function(c) rowSums(Q*c[cells,,drop=FALSE]))
  pred <- lapply(names(raw),function(j) {
    bounds <- if (j=="pi") fit$schedule$propensity_bounds else fit$schedule$outcome_bounds
    pmin(pmax(raw[[j]],bounds[1]),bounds[2])
  }) |> setNames(names(raw))
  pred$raw <- raw
  pred$clip_fraction <- vapply(names(raw),function(j) mean(raw[[j]]!=pred[[j]]),numeric(1)) |> setNames(names(raw))
  pred
}

affine_crossfit <- function(design,T,Y,seed,folds=5L,full_stats=NULL,timing=FALSE) {
  if (folds!=5L || design$n<folds) stop("Five nonempty folds are required")
  fold <- affine_with_seed(seed,sample(rep(seq_len(folds),length.out=design$n)))
  if(is.null(full_stats)) full_stats <- affine_stats(design,T,Y)
  fits <- lapply(seq_len(folds),function(k) {
    removed <- affine_stats(design,T,Y,which(fold==k))
    tryCatch(affine_fit_from_stats(affine_stats_difference(full_stats,removed),design$schedule,which(fold!=k),timing),
      affine_component_error=function(e) {e$fold <- as.integer(k); stop(e)})
  })
  keys <- names(full_stats); out <- setNames(rep(list(numeric(design$n)),length(keys)),keys)
  raw <- out
  # Row-only prediction avoids storing repeated n-by-four fit prediction arrays.
  for (k in seq_len(folds)) {
    ix <- which(fold==k)
    for (j in keys) {
      val <- rowSums(design$Q[ix,,drop=FALSE]*fits[[k]]$coefficients[[j]][design$cell[ix],,drop=FALSE])
      raw[[j]][ix] <- val
      bounds <- if(j=="pi") design$schedule$propensity_bounds else design$schedule$outcome_bounds
      out[[j]][ix] <- pmin(pmax(val,bounds[1]),bounds[2])
    }
  }
  out$raw <- raw
  out$clip_fraction <- vapply(keys,function(j) mean(raw[[j]]!=out[[j]]),numeric(1)) |> setNames(keys)
  list(folds=fold,fits=fits,predictions=out)
}

affine_delete <- function(fit,design,T,Y,indices) {
  indices <- affine_validate_retained(indices,design$n)
  if (length(indices)>2 || any(!indices%in%fit$retain)) stop("Deletion must select at most two retained rows")
  sch <- fit$schedule; penalty <- diag(sch$n*sch$lambda/sch$J,sch$q_dim)
  for (i in indices) {
    cell <- design$cell[i]; q <- design$Q[i,]
    for (j in names(fit$coefficients)) {
      D <- switch(j,pi=1,m=1,mu0=1-T[i],mu1=T[i])
      if (D==0) next
      R <- if(j=="pi") T[i] else Y[i]
      st <- fit$stats[[j]]; A <- st$S[,,cell]+penalty
      v <- affine_solve_block(A,q); denom <- 1-sum(q*v)
      st$S[,,cell] <- st$S[,,cell]-tcrossprod(q)
      st$h[,cell] <- st$h[,cell]-q*R; st$count[cell] <- st$count[cell]-1L
      if (!is.finite(denom) || denom<1e-8) {
        fit$coefficients[[j]][cell,] <- affine_solve_block(st$S[,,cell]+penalty,st$h[,cell])
        fit$update_fallbacks <- fit$update_fallbacks+1L
      } else {
        c <- fit$coefficients[[j]][cell,]
        fit$coefficients[[j]][cell,] <- c-v*(R-sum(q*c))/denom
      }
      fit$stats[[j]] <- st
    }
  }
  fit$retain <- setdiff(fit$retain,indices)
  fit$diagnostics <- affine_gram_diagnostics(fit$stats,sch)
  fit$coefficient_norms <- vapply(fit$coefficients,function(c) sqrt(sum(c^2)/sch$J),numeric(1))
  fit
}

affine_population_fit <- function(schedule,design_id=schedule$design_id,ridge=FALSE) {
  list(coefficients=affine_population_coefficients(schedule,design_id,ridge),schedule=schedule,
       parameterization="cell_q",method=if(ridge) "population_ridge" else "population_projection")
}
