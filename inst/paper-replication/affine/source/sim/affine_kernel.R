# Product Gaussian localization on the reporting coordinates only. Bandwidths
# use the specified reporting population, never the realized draw or coverage.
affine_kernel_population <- function(n,h,x=c(.5,.5),design_id="C0") {
  if(length(h)!=1L || !is.finite(h) || h<=0 || length(x)!=2L ||
     any(!is.finite(x)) || any(x<=0 | x>=1)) stop("Invalid kernel population geometry")
  sparse <- affine_model_profile(design_id)$sparse
  axis_mass <- function(power) {
    scale <- sqrt(power)/h
    whole <- pnorm((1-x)*scale)-pnorm(-x*scale)
    fraction <- if(sparse) 1.5*whole-(pnorm((.75-x)*scale)-pnorm((.25-x)*scale)) else whole
    sqrt(2*pi)*h/sqrt(power)*fraction
  }
  # Keep the original uniform expressions bitwise unchanged for old runs.
  mass <- if(!sparse) prod(sqrt(2*pi)*h*(pnorm((1-x)/h)-pnorm(-x/h))) else prod(axis_mass(1))
  square_mass <- if(!sparse) prod(sqrt(pi)*h*(pnorm(sqrt(2)*(1-x)/h)-pnorm(-sqrt(2)*x/h))) else prod(axis_mass(2))
  list(mass=mass,square_mass=square_mass,N=n*mass^2/square_mass)
}

affine_kernel_bandwidth <- function(n,s,design_id="C0") {
  reference <- 1/sum(affine_rank_weights(n,s)^2)
  if(s<=1 || reference>=n) stop("Population-matched local kernel requires s>1")
  initial <- sqrt(reference/(4*pi*n))
  gap <- function(h) affine_kernel_population(n,h,design_id=design_id)$N/reference-1
  upper <- initial
  while(gap(upper)<0) upper <- 2*upper
  uniroot(gap,c(initial/2,upper),tol=1e-12)$root
}

affine_kernel_weights <- function(X,x=c(.5,.5),h,s,design_id="C0") {
  affine_unit_matrix(X,2,"X"); n <- nrow(X)
  if(n<2 || length(x)!=2 || any(!is.finite(x)) || any(x<=0 | x>=1) ||
     length(h)!=1 || !is.finite(h) || h<=0) stop("Invalid kernel bandwidth/reporting point")
  distance <- sqrt(rowSums(sweep(X,2,x)^2))
  log_kernel <- -.5*(distance/h)^2
  if(any(!is.finite(log_kernel))) stop("Kernel exponent overflow")
  kernel <- exp(log_kernel-max(log_kernel)); weights <- kernel/sum(kernel)
  ord <- order(distance,seq_len(n))
  radial <- function(p) distance[ord[which(cumsum(weights[ord])>=p)[1L]]]
  zero_count <- sum(kernel==0)
  # Conservative bound for any underflowed mass after the common log shift.
  log_bound <- if(zero_count) log(zero_count)+log(.Machine$double.xmin)-log(sum(kernel)) else -Inf
  list(method="kernel_gaussian",n=n,s=s,x=x,h=h,weights=weights,distance=distance,
       N=1/sum(weights^2),N_population=affine_kernel_population(n,h,x,design_id)$N,
       mathematical_support=n,numerical_support=sum(weights>0),
       omitted_mass_bound=exp(log_bound),log_omitted_mass_bound=log_bound,
       normalization_error=abs(sum(weights)-1),weighted_radius=sum(weights*distance),
       radius90=radial(.9),radius95=radial(.95),largest_weight=max(weights))
}

affine_kernel_joint_pseudo <- function(equations,local) {
  PHI <- cbind(phi=equations$phi,B=equations$slope$B,b=equations$slope$b)
  w <- local$weights; retained <- 1-w
  dominant <- which(retained<sqrt(.Machine$double.eps))
  base <- drop(crossprod(w,PHI))
  # Same h, nuisances and fitted target; only remove a weight and renormalize.
  deleted <- (matrix(base,local$n,3L,byrow=TRUE)-PHI*w)/pmax(retained,sqrt(.Machine$double.eps))
  for(j in dominant) {
    # Avoid cancellation in both the retained mass and retained numerator.
    kept <- w; kept[j] <- 0; mass <- sum(kept)
    if(mass<=0) stop("No numerically retained kernel mass after deletion")
    deleted[j,] <- drop(crossprod(kept/mass,PHI))
  }
  Q <- local$n*matrix(base,local$n,3L,byrow=TRUE)-(local$n-1L)*deleted
  list(base=base,deleted=deleted,Q=Q,n=local$n,deletion_recomputations=length(dominant))
}

affine_kernel_population_localizer <- function(h,design_id="C0") {
  affine_kernel_population(1,h,design_id=design_id)
  profile <- affine_model_profile(design_id)
  if(profile$base!="C2") return(list(target=1/8,drift=0,integration_error=0,method="kernel_reflection_symmetry"))
  if(profile$sparse) {
    end <- .5/h; cap <- min(end,12); split <- min(.25/h,cap)
    cuts <- unique(c(0,split,cap))
    values <- lapply(seq_len(length(cuts)-1L),function(i) integrate(function(t)
      affine_reporting_density(.5+h*t,design_id)*exp(-t^2/2)*2*sin(4*pi*h*t)^2,
      cuts[i],cuts[i+1],rel.tol=1e-10,abs.tol=1e-12))
    denominator <- sqrt(affine_kernel_population(1,h,design_id=design_id)$mass)/(2*h)
    tail <- if(cap<end) 3*sqrt(2*pi)*(pnorm(cap,lower.tail=FALSE)-pnorm(end,lower.tail=FALSE)) else 0
    drift <- sum(vapply(values,`[[`,numeric(1),"value"))/(8*denominator)
    error <- (sum(vapply(values,`[[`,numeric(1),"abs.error"))+tail)/(8*denominator)
    return(list(target=1/8+drift,drift=drift,integration_error=error,
      omitted_tail_bound=tail/(8*denominator),method="kernel_sparse_density_quadrature"))
  }
  # Standardized integration avoids a narrow peak at small bandwidths. The
  # exact normalizer includes the square boundary; omitted tail is bounded.
  end <- .5/h; cap <- min(end,12)
  denominator <- sqrt(2*pi)*(pnorm(end)-.5)
  integral <- integrate(function(t) exp(-t^2/2)*2*sin(4*pi*h*t)^2,
        0,cap,rel.tol=1e-10,abs.tol=1e-12)
  tail <- if(cap<end) 2*sqrt(2*pi)*(pnorm(cap,lower.tail=FALSE)-pnorm(end,lower.tail=FALSE)) else 0
  drift <- integral$value/(8*denominator)
  list(target=1/8+drift,drift=drift,
       integration_error=(integral$abs.error+tail)/(8*denominator),
       omitted_tail_bound=tail/(8*denominator),method="kernel_truncated_gaussian_quadrature")
}
