# Fresh laws condition on X only; Z,T and outcome noise remain integrated out.
.affine_integral_cache <- new.env(parent=emptyenv())

affine_sine_antiderivative <- function(t,k) {
  if(k==0) return(t)
  if(k==1) return(-cos(t))
  -sin(t)^(k-1)*cos(t)/k+(k-1)/k*affine_sine_antiderivative(t,k-2)
}

affine_integrals <- function(M) {
  key <- as.character(M)
  if(exists(key,.affine_integral_cache,inherits=FALSE)) return(get(key,.affine_integral_cache))
  lo <- (0:(M-1))/M; hi <- (1:M)/M
  moments <- matrix(vapply(0:5,function(k) M*(affine_sine_antiderivative(2*pi*hi,k)-affine_sine_antiderivative(2*pi*lo,k))/(2*pi),numeric(M)),M,6L)
  reciprocal <- lapply(seq_len(M),function(i) integrate(function(z)
       1/(.5+sin(2*pi*z)/8)+1/(.5-sin(2*pi*z)/8),lo[i],hi[i],rel.tol=1e-11,abs.tol=1e-12))
  cells <- affine_cell_moments(M)
  fourth <- cells$bins[,4]
  out <- c(cells,list(s4_moments=moments[fourth,,drop=FALSE],
       reciprocal=vapply(reciprocal,`[[`,numeric(1),"value")[fourth]*M,
       integration_error=max(vapply(reciprocal,`[[`,numeric(1),"abs.error"))*M))
  assign(key,out,.affine_integral_cache); out
}

affine_poly_add <- function(A,B) {
  out <- matrix(0,nrow(A),max(ncol(A),ncol(B)))
  out[,seq_len(ncol(A))] <- A
  out[,seq_len(ncol(B))] <- out[,seq_len(ncol(B)),drop=FALSE]+B
  out
}

affine_poly_mul <- function(A,B) {
  out <- matrix(0,nrow(A),ncol(A)+ncol(B)-1L)
  for(i in seq_len(ncol(A))) for(j in seq_len(ncol(B))) out[,i+j-1L] <- out[,i+j-1L]+A[,i]*B[,j]
  out
}

affine_fresh <- function(fit,u,design_id=fit$schedule$design_id,a0=1/8,components=TRUE) {
  if(length(u)!=2 || any(!is.finite(u)) || any(u<0 | u>1)) stop("Invalid fresh reporting point")
  if(identical(fit$method,"oracle")) return(list(mse=setNames(numeric(4),c("pi","mu0","mu1","m")),
      cate=list(signal_bias=0,score_mean=0),slope=list(coefficient_increment=0,product=0,local_remainder=0,mixed_spatial=0,score_mean=0),
      component_mse=c(cate_b=0,slope_B=0,slope_b=0,slope_score=0),integration_error=0))
  C <- affine_integrals(fit$schedule$M); S <- C$sine_mean; V <- C$sine_variance
  J <- fit$schedule$J; h <- if(affine_model_profile(design_id)$base=="C1") 1/8 else 0
  noise <- affine_noise_cell_variance(u,C$bins,fit$schedule$M,design_id)
  Z <- (C$bins-.5)/fit$schedule$M; X <- matrix(u,J,2,byrow=TRUE)
  pred <- affine_predict(fit,X,Z); e <- pred$pi
  theta <- affine_truth(matrix(u,1,2),matrix(0,1,4),design_id)$theta
  gp <- (S[,1]+u[1]*S[,2]+u[2]*S[,3])/8
  vg <- (V[,1]+u[1]^2*V[,2]+u[2]^2*V[,3])/64
  barp <- .5+S[,4]/8; varp <- V[,4]/64
  barg <- gp+h*S[,4]; barm <- barg+barp*theta
  truth_mean <- list(pi=barp,mu0=barg,mu1=barg+theta,m=barm)
  truth_var <- list(pi=varp,mu0=vg+h^2*V[,4],mu1=vg+h^2*V[,4],m=vg+(h+theta/8)^2*V[,4])
  mse <- vapply(names(truth_mean),function(j) mean((pred[[j]]-truth_mean[[j]])^2+truth_var[[j]]),numeric(1))
  covpg <- h*V[,4]/8
  product <- mean((e-barp)*(pred$m-barm)+covpg+theta*varp)
  signal_bias <- mean(((e-barp)*(pred$mu1-barg-theta)+covpg)/e+
                       ((e-barp)*(pred$mu0-barg)+covpg)/(1-e))
  local_remainder <- theta*mse["pi"]-product
  mixed <- (a0-theta)*mse["pi"]
  out <- list(mse=mse,cate=list(signal_bias=signal_bias,score_mean=-signal_bias),
    slope=list(coefficient_increment=unname(mse["pi"]),product=product,local_remainder=unname(local_remainder),
               mixed_spatial=unname(mixed),score_mean=unname(local_remainder+mixed)))
  if(!components) return(out)
  integrate_poly <- function(A) rowSums(A*C$s4_moments[,seq_len(ncol(A)),drop=FALSE])
  norms <- c(cate_b=0,slope_B=0,slope_b=0,slope_score=0)
  ppoly <- cbind(rep(.5,J),rep(1/8,J))
  for(t in 0:1) {
    prob <- if(t==1) ppoly else cbind(rep(.5,J),rep(-1/8,J))
    v <- t-e; v0 <- cbind(rep(t-.5,J),rep(-1/8,J))
    v0sq <- affine_poly_mul(v0,v0)
    dB <- affine_poly_add(matrix(v^2,J,1),-v0sq)
    # b_fit-b_oracle: its conditional mean is polynomial in sin(2*pi Z4).
    db <- affine_poly_add(cbind(-v*(gp+t*theta-pred$m),-v*h),theta*v0sq)
    de <- cbind(e-.5,rep(-1/8,J))
    conditional_var <- affine_poly_add(matrix(v^2*vg,J,1),affine_poly_mul(de,de)*noise)
    ds <- affine_poly_add(a0*dB,db)
    norms["slope_B"] <- norms["slope_B"]+mean(integrate_poly(affine_poly_mul(prob,affine_poly_mul(dB,dB))))
    norms["slope_b"] <- norms["slope_b"]+mean(integrate_poly(affine_poly_mul(prob,affine_poly_add(affine_poly_mul(db,db),conditional_var))))
    norms["slope_score"] <- norms["slope_score"]+mean(integrate_poly(affine_poly_mul(prob,affine_poly_add(affine_poly_mul(ds,ds),conditional_var))))
    cf <- if(t==1) 1/e else -1/(1-e)
    mu <- if(t==1) pred$mu1 else pred$mu0
    dphi <- cbind(pred$mu1-pred$mu0-theta+cf*(gp+t*theta-mu),cf*h)
    moment <- affine_poly_add(affine_poly_mul(dphi,dphi),matrix(cf^2*vg,J,1))
    norms["cate_b"] <- norms["cate_b"]+mean(integrate_poly(affine_poly_mul(prob,moment)))
  }
  # Integrate inverse propensity terms separately with recorded quadrature error.
  norms["cate_b"] <- norms["cate_b"]+mean(noise*(barp/e^2+(1-barp)/(1-e)^2-2/e-2/(1-e)+C$reciprocal))
  norms <- pmax(norms,0)
  out$component_mse <- norms; out$integration_error <- C$integration_error*max(noise)
  out
}

affine_mean_grid <- function(fit,U,design_id=fit$schedule$design_id,a0=1/8,chunk=512L) {
  affine_unit_matrix(U,2,"U")
  mean_keys <- c("cate","slope","slope_coefficient_increment","slope_product","slope_local_remainder","slope_mixed_spatial")
  if(identical(fit$method,"oracle")) return(setNames(rep(list(numeric(nrow(U))),length(mean_keys)),mean_keys))
  C <- affine_cell_moments(fit$schedule$M); S <- C$sine_mean; V <- C$sine_variance
  J <- fit$schedule$J; h <- if(affine_model_profile(design_id)$base=="C1") 1/8 else 0
  out <- setNames(rep(list(numeric(nrow(U))),length(mean_keys)),mean_keys)
  for(start in seq.int(1L,nrow(U),by=chunk)) {
    ix <- start:min(nrow(U),start+chunk-1L); u <- U[ix,,drop=FALSE]; nr <- nrow(u)
    Q <- affine_basis(u,fit$schedule$q_dim)
    pp <- lapply(names(fit$coefficients),function(j) {
      val <- Q%*%t(fit$coefficients[[j]]); bounds <- if(j=="pi") c(.25,.75) else c(-1,1)
      pmin(pmax(val,bounds[1]),bounds[2])
    }) |> setNames(names(fit$coefficients))
    repcol <- function(v) matrix(v,nr,J,byrow=TRUE)
    theta <- affine_truth(u,matrix(0,nr,4),design_id)$theta
    theta <- matrix(theta,nr,J)
    barg <- repcol(S[,1]/8+h*S[,4])+tcrossprod(u[,1],S[,2]/8)+tcrossprod(u[,2],S[,3]/8)
    barp <- repcol(.5+S[,4]/8); varp <- repcol(V[,4]/64); covpg <- repcol(h*V[,4]/8)
    de <- pp$pi-barp
    mse <- rowMeans(de^2+varp)
    product <- rowMeans(de*(pp$m-barg-barp*theta)+covpg+theta*varp)
    bias <- rowMeans((de*(pp$mu1-barg-theta)+covpg)/pp$pi+(de*(pp$mu0-barg)+covpg)/(1-pp$pi))
    out$cate[ix] <- -bias; out$slope[ix] <- a0*mse-product
    out$slope_coefficient_increment[ix] <- mse; out$slope_product[ix] <- product
    out$slope_local_remainder[ix] <- theta[,1]*mse-product
    out$slope_mixed_spatial[ix] <- (a0-theta[,1])*mse
  }
  out
}

affine_covariance_transfer <- function(fitted,oracle,application,target=1/8,N) {
  if(is.null(fitted$joint) || is.null(oracle$joint)) return(list(
       transfer=NA_real_,component=NA_real_,target=NA_real_,cross=NA_real_,
       point_oracle_scaled=NA_real_,oracle_root_scaled=NA_real_,
       oracle_Gamma=oracle[[application]]$Gamma_hat,status="unavailable"))
  qa <- affine_frozen_pseudo(fitted$joint,fitted[[application]]$estimate,application)
  qo <- affine_frozen_pseudo(oracle$joint,target,application)
  component <- affine_frozen_pseudo(fitted$joint,target,application)-qo
  movement <- qa-affine_frozen_pseudo(fitted$joint,target,application)
  center <- function(v) v-mean(v)
  factor <- N/(length(qa)*(length(qa)-1))
  component <- center(component); movement <- center(movement)
  list(transfer=factor*sum((component+movement)^2),component=factor*sum(component^2),
       target=factor*sum(movement^2),cross=2*factor*sum(component*movement),
       point_oracle_scaled=factor*sum(center(qo)^2),oracle_root_scaled=oracle[[application]]$score_cov_scaled["jackknife"],
       oracle_Gamma=oracle[[application]]$Gamma_hat)
}

affine_population_localizer <- function(s,design_id="C0") {
  if(s<1 || s!=floor(s)) stop("Invalid localizer reference")
  profile <- affine_model_profile(design_id)
  if(profile$sparse) return(affine_sparse_population_localizer(s))
  if(profile$base!="C2") return(list(target=1/8,drift=0,integration_error=0,method="reflection_symmetry"))
  if(s==1) return(list(target=1/4,drift=1/8,integration_error=0,method="uniform_integral"))
  key <- paste0("localizer:C2:",s)
  if(exists(key,.affine_integral_cache,inherits=FALSE)) return(get(key,.affine_integral_cache))
  A <- 1/4; cap <- min(40,s*pi/4)
  # t=s*pi*r^2 concentrates the integral at its shrinking interior radius.
  central <- integrate(function(t) (1-t/s)^(s-1)*A/2*(1-besselJ(8*pi*sqrt(t/(s*pi)),0)),
        0,cap,rel.tol=1e-10,abs.tol=1e-12)
  tail <- if(cap<s*pi/4) A*(1-cap/s)^s else A*(1-pi/4)^s
  outside <- 0; error <- central$abs.error
  if(cap==s*pi/4 && tail>1e-12) {
    angular_error <- 0
    integrand <- function(radii) vapply(radii,function(r) {
      lower <- acos(.5/r); upper <- pi/2-lower
      angle <- integrate(function(a) A/4*(2-cos(8*pi*r*cos(a))-cos(8*pi*r*sin(a))),
            lower,upper,rel.tol=1e-10,abs.tol=1e-13)
      angular_error <<- max(angular_error,angle$abs.error)
      F <- pi*r^2-4*(r^2*acos(.5/r)-.5*sqrt(r^2-.25))
      4*s*max(0,1-F)^(s-1)*r*angle$value
    },numeric(1))
    outer <- integrate(integrand,.5,sqrt(.5),rel.tol=1e-10,abs.tol=1e-12)
    outside <- outer$value
    error <- error+outer$abs.error+4*s*sqrt(.5)*(sqrt(.5)-.5)*angular_error
    tail <- 0
  }
  drift <- central$value+outside
  out <- list(target=1/8+drift,drift=drift,integration_error=error+tail,
              omitted_tail_bound=tail,method="radial_bessel_and_square_boundary_quadrature")
  assign(key,out,.affine_integral_cache); out
}

affine_sparse_population_localizer <- function(s) {
  key <- paste0("localizer:C2S:",s)
  if(exists(key,.affine_integral_cache,inherits=FALSE)) return(get(key,.affine_integral_cache))
  top <- s*pi/64; cap <- min(40,top)
  central <- integrate(function(t) (1-t/s)^(s-1)/8*(1-besselJ(8*pi*sqrt(t/(s*pi/4)),0)),
    0,cap,rel.tol=1e-10,abs.tol=1e-12)
  tail <- (1-cap/s)^s/4; outside <- 0; error <- central$abs.error
  if(cap==top && tail>1e-12) {
    angle_error <- 0
    integrand <- function(radii) vapply(radii,function(r) {
      lo <- if(r>.5) acos(.5/r) else 0; hi <- pi/2-lo
      cuts <- sort(unique(c(lo,hi,asin(min(1,.25/r)),acos(min(1,.25/r)))))
      cuts <- cuts[cuts>=lo & cuts<=hi]
      angles <- lapply(seq_len(length(cuts)-1L),function(i) integrate(function(a)
        affine_reporting_density(.5+r*cos(a),"C2S")*affine_reporting_density(.5+r*sin(a),"C2S")*
          (sin(4*pi*r*cos(a))^2+sin(4*pi*r*sin(a))^2)/8,
        cuts[i],cuts[i+1],rel.tol=1e-10,abs.tol=1e-13))
      angle_error <<- max(angle_error,sum(vapply(angles,`[[`,numeric(1),"abs.error")))
      4*s*(1-affine_reporting_ball(r,"C2S"))^(s-1)*r*sum(vapply(angles,`[[`,numeric(1),"value"))
    },numeric(1))
    cuts <- c(.25,sqrt(2)/4,.5,sqrt(.5))
    outer <- lapply(1:3,function(i) integrate(integrand,cuts[i],cuts[i+1],rel.tol=1e-9,abs.tol=1e-12))
    outside <- sum(vapply(outer,`[[`,numeric(1),"value"))
    error <- error+sum(vapply(outer,`[[`,numeric(1),"abs.error"))+s*angle_error
    tail <- 0
  }
  drift <- central$value+outside
  out <- list(target=1/8+drift,drift=drift,integration_error=error+tail,
    omitted_tail_bound=tail,method="sparse_radial_and_piecewise_angular_quadrature")
  assign(key,out,.affine_integral_cache); out
}

affine_selected_deletions <- function(data,fit,local,singletons=16L,pairs=8L) {
  design <- affine_design(data$X,data$Z,fit$schedule); n <- data$n
  leverage <- numeric(n)
  penalty <- diag(fit$schedule$n*fit$schedule$lambda/fit$schedule$J,fit$schedule$q_dim)
  for(c in unique(design$cell)) {
    ix <- which(design$cell==c); Q <- design$Q[ix,,drop=FALSE]
    for(j in names(fit$coefficients)) {
      v <- rowSums(Q*(Q%*%solve(fit$stats[[j]]$S[,,c]+penalty)))
      mask <- switch(j,pi=rep(1,length(ix)),m=rep(1,length(ix)),mu0=1-data$T[ix],mu1=data$T[ix])
      leverage[ix] <- pmax(leverage[ix],v*mask)
    }
  }
  reps <- unlist(lapply(split(seq_len(n),interaction(design$cell,data$T,drop=TRUE)),head,1L))
  candidate_groups <- list(head(seq_len(n),4),head(order(-local$weights),4),head(order(-leverage),4),head(reps,4))
  candidates <- unlist(candidate_groups)
  reasons <- rep(c("fixed_id","large_local_weight","large_leverage","cell_arm_representative"),lengths(candidate_groups))
  keep <- !duplicated(candidates); indices <- head(candidates[keep],singletons)
  log <- data.frame(index=indices,reason=head(reasons[keep],length(indices)))
  combos <- if(length(indices)>1) t(combn(indices,2)) else matrix(integer(),0,2)
  if(nrow(combos)) {
    same <- design$cell[combos[,1]]==design$cell[combos[,2]]
    preferred <- unique(c(which(same)[1],which(!same)[1],seq_len(nrow(combos))))
    preferred <- preferred[!is.na(preferred)]
    combos <- combos[head(preferred,pairs),,drop=FALSE]
  }
  rL <- if(identical(local$method,"kernel_gaussian")) local$h*sqrt(8*log(n)) else sqrt(8*log(n)/local$s)
  r <- min(.25,rL/2)
  grid <- rbind(local$x,local$x+c(r,0),local$x-c(r,0),local$x+c(0,r),local$x-c(0,r))
  C <- affine_cell_moments(fit$schedule$M); Z <- (C$bins-.5)/fit$schedule$M
  Xeval <- grid[rep(seq_len(5),each=nrow(Z)),,drop=FALSE]
  Zeval <- Z[rep(seq_len(nrow(Z)),5),,drop=FALSE]
  base <- affine_predict(fit,Xeval,Zeval)
  grid_errors <- function(f) t(vapply(seq_len(nrow(grid)),function(i)
    affine_fresh(f,grid[i,],data$design_id,components=FALSE)$mse,numeric(4)))
  maximum_errors <- grid_errors(fit)
  widths <- c(raw=0,projected=0,cate_b=0,slope_B=0,slope_b=0,slope_score=0)
  fallback <- 0L
  selected <- c(lapply(indices,function(i) i),lapply(seq_len(nrow(combos)),function(i) combos[i,]))
  for(ix in selected) {
    deleted <- affine_delete(fit,design,data$T,data$Y,ix); fallback <- fallback+deleted$update_fallbacks
    maximum_errors <- pmax(maximum_errors,grid_errors(deleted))
    p <- affine_predict(deleted,Xeval,Zeval)
    widths["raw"] <- max(widths["raw"],max(abs(unlist(p$raw)-unlist(base$raw))))
    keys <- c("pi","mu0","mu1","m")
    widths["projected"] <- max(widths["projected"],max(abs(unlist(p[keys])-unlist(base[keys]))))
    for(t in 0:1) for(y in c(-1,1)) {
      W <- list(n=nrow(Xeval),T=rep(t,nrow(Xeval)),Y=rep(y,nrow(Xeval)))
      ep <- affine_equations(W,p); eb <- affine_equations(W,base)
      db <- ep$slope$b-eb$slope$b; dB <- ep$slope$B-eb$slope$B
      current <- c(cate_b=max(abs(ep$cate$b-eb$cate$b)),slope_B=max(abs(dB)),slope_b=max(abs(db)),slope_score=max(abs(dB/8+db)))
      widths[names(current)] <- pmax(widths[names(current)],current)
    }
  }
  list(singletons=indices,selection_log=log,pairs=combos,grid=grid,rL=rL,
       grid_maximum_mse=maximum_errors,grid_prediction_radius_squared=max(rowSums(maximum_errors)),
       widths=widths,scaled_widths=sqrt(local$N)*widths,update_fallbacks=fallback,
       order_reference_K_squared_over_n=fit$schedule$K^2/fit$schedule$n,
       scaled_order_reference=sqrt(local$N)*fit$schedule$K^2/fit$schedule$n,
       order_reference_scope="unspecified_constant_order_reference_not_a_numerical_bound",
       scope="selected_states_and_five_point_grid_not_proof_supremum")
}

affine_diagnostics <- function(data,fit,equations,local,design_id,diagnostic_spec=list(),inference=NULL) {
  a0 <- 1/8; intensive <- isTRUE(diagnostic_spec$intensive)
  fitted <- if(is.null(inference)) affine_inference(equations,local) else inference
  oracle_eq <- affine_equations(data,data$truth); oracle <- affine_inference(oracle_eq,local)
  keys <- c("pi","mu0","mu1","m")
  if(!is.null(fit$fits)) {
    fresh <- lapply(fit$fits,affine_fresh,u=local$x,design_id=design_id)
    masses <- vapply(seq_along(fit$fits),function(k) sum(local$weights[fit$folds==k]),numeric(1))
    accuracy <- list(fold_average=Reduce(`+`,lapply(fresh,`[[`,"mse"))/5,
      local_mass_average=Reduce(`+`,Map(function(f,w) f$mse*w,fresh,masses)),fold_masses=masses,fresh_functions=fresh)
    predictions <- fit$predictions
  } else {
    accuracy <- affine_fresh(fit,local$x,design_id)
    predictions <- if(identical(fit$method,"oracle")) data$truth else affine_predict(fit,data$X,data$Z)
  }
  training_mse <- vapply(keys,function(j) mean((predictions[[j]]-data$truth[[j]])^2),numeric(1))
  v0 <- data$T-data$truth$pi; Gamma0 <- sum(local$weights*v0^2)
  weighted_truth <- c(cate=sum(local$weights*data$truth$theta),slope=sum(local$weights*v0^2*data$truth$theta)/Gamma0)
  covariance <- lapply(c("cate","slope"),function(app) affine_covariance_transfer(fitted,oracle,app,a0,local$N)) |> setNames(c("cate","slope"))
  out <- list(accuracy=accuracy,training_mse=training_mse,weighted_truth=weighted_truth,
      population_localizer=if(identical(local$method,"kernel_gaussian"))
        affine_kernel_population_localizer(local$h,design_id) else affine_population_localizer(local$s,design_id),
      spatial_drift=weighted_truth-a0,covariance=covariance,
      oracle_constants=c(Gamma=31/128,Sigma_cate=1/(4*sqrt(15)),Sigma_slope=31/8192,Omega_slope=2/31),
      conditional_noise=if(!affine_model_profile(design_id)$heteroskedastic)
        c(cate_X=1/(4*sqrt(15)*local$N),
          cate_VT=sum(local$weights^2*(data$T/data$truth$pi^2+(1-data$T)/(1-data$truth$pi)^2))/64,
          slope_VT=sum(local$weights^2*v0^2)/(64*Gamma0^2)) else
        c(cate_X=16/sqrt(15)*sum(local$weights^2*affine_conditional_noise_variance(data$X[,1],design_id)),
          cate_VT=sum(local$weights^2*(data$T/data$truth$pi^2+(1-data$T)/(1-data$truth$pi)^2)*data$truth$noise_scale^2),
          slope_VT=sum(local$weights^2*v0^2*data$truth$noise_scale^2)/Gamma0^2))
  out$clip_frequency <- if(is.null(predictions$clip_fraction)) setNames(numeric(4),keys) else predictions$clip_fraction
  out$kernel_deletion_recomputations <- if(identical(local$method,"kernel_gaussian") && !is.null(fitted$joint))
    fitted$joint$deletion_recomputations else NA_integer_
  out$gram <- if(!is.null(fit$fits)) lapply(fit$fits,`[[`,"diagnostics") else fit$diagnostics
  out$local <- local[c("method","h","N","N_population","s","weighted_radius","radius90","radius95","largest_weight",
     "mathematical_support","numerical_support","omitted_mass_bound","log_omitted_mass_bound","normalization_error")]
  if(intensive) {
    if(!is.null(fit$fits)) {
      mean_keys <- c("cate","slope","slope_coefficient_increment","slope_product","slope_local_remainder","slope_mixed_spatial")
      means <- setNames(rep(list(numeric(data$n)),length(mean_keys)),mean_keys)
      for(k in seq_along(fit$fits)) {
        ix <- which(fit$folds==k); fm <- affine_mean_grid(fit$fits[[k]],data$X[ix,,drop=FALSE],design_id)
        for(app in names(means)) means[[app]][ix] <- fm[[app]]
      }
    } else means <- affine_mean_grid(fit,data$X,design_id)
    decomposition <- lapply(c("cate","slope"),function(app) {
      delta <- (equations[[app]]$B-oracle_eq[[app]]$B)*a0+equations[[app]]$b-oracle_eq[[app]]$b
      meanpart <- sum(local$weights*means[[app]])
      centered <- sum(local$weights*(delta-means[[app]]))
      gap <- fitted[[app]]$estimate-oracle[[app]]$estimate
      spatial <- if(app=="slope") (fitted$slope$Gamma_hat-oracle$slope$Gamma_hat)*(oracle$slope$estimate-a0) else 0
      identity <- fitted[[app]]$Gamma_hat*gap+meanpart+centered+spatial
      out_parts <- list(mean_part=meanpart,centered_part=centered,coefficient_spatial=spatial,
           identity_error=identity,scaled_mean=sqrt(local$N)*meanpart,scaled_centered=sqrt(local$N)*centered)
      for(key in c("coefficient_increment","product","local_remainder","mixed_spatial"))
        out_parts[[paste0("weighted_",key)]] <- if(app=="slope") sum(local$weights*means[[paste0("slope_",key)]]) else NA_real_
      out_parts
    }) |> setNames(c("cate","slope"))
    out$decomposition <- decomposition
    if(identical(fit$method,"sample_ridge")) out$deletions <- affine_selected_deletions(data,fit,local,
       if(is.null(diagnostic_spec$singletons)) 16L else diagnostic_spec$singletons,
       if(is.null(diagnostic_spec$pairs)) 8L else diagnostic_spec$pairs)
  }
  out
}
