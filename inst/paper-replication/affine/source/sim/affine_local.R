# k=1 exact DNN weights and one joint frozen-array covariance engine call.
affine_rank_weights <- function(n,s) {
  if (length(n)!=1 || length(s)!=1 || n<1 || s<1 || s>n || n!=floor(n) || s!=floor(s)) stop("Invalid rank law")
  w <- numeric(n); w[1] <- s/n; support <- n-s+1L
  if (support>1L) for(r in seq_len(support-1L)) w[r+1L] <- w[r]*(n-r-s+1)/(n-r)
  w
}

affine_weights <- function(X,x=c(.5,.5),s) {
  affine_unit_matrix(X,2,"X"); n <- nrow(X)
  if (length(x)!=2 || any(!is.finite(x)) || any(x<0 | x>1) || s>=n || s<1 || s!=floor(s)) stop("Invalid reporting point/order; same s must fit n-1")
  distance <- sqrt(rowSums(sweep(X,2,x)^2))
  ord <- order(distance,seq_len(n)); ranks <- integer(n); ranks[ord] <- seq_len(n)
  wn <- affine_rank_weights(n,s); wnm1 <- affine_rank_weights(n-1L,s)
  weights <- wn[ranks]; support <- n-s+1L
  zeros <- which(wn[seq_len(support)]==0)
  logmass <- if (length(zeros)) lchoose(n-zeros[1]+1,s)-lchoose(n,s) else -Inf
  radial <- function(p) distance[ord[which(cumsum(wn)>=p)[1L]]]
  list(method="dnn",n=n,s=s,x=x,weights=weights,ranks=ranks,order=ord,distance=distance,
       rank_weights=wn,deleted_rank_weights=wnm1,N=1/sum(wn^2),
       mathematical_support=support,numerical_support=sum(wn>0),
       omitted_mass_bound=exp(logmass),log_omitted_mass_bound=logmass,
       normalization_error=abs(sum(wn)-1),weighted_radius=sum(weights*distance),
       radius90=radial(.9),radius95=radial(.95),largest_weight=max(weights))
}

affine_equations <- function(data,predictions) {
  n <- data$n; keys <- c("pi","mu0","mu1","m")
  if (!all(keys%in%names(predictions)) || any(vapply(predictions[keys],function(x)
      length(x)!=n || any(!is.finite(x)),logical(1)))) stop("Missing or nonfinite nuisance predictions")
  e <- predictions$pi
  if (any(e<.25 | e>.75)) stop("Score propensity must use fixed overlap projection")
  T <- data$T; Y <- data$Y
  phi <- predictions$mu1-predictions$mu0+T*(Y-predictions$mu1)/e-(1-T)*(Y-predictions$mu0)/(1-e)
  v <- T-e
  list(n=n,phi=phi,cate=list(B=rep(1,n),b=-phi),slope=list(B=v^2,b=-v*(Y-predictions$m)))
}

affine_joint_pseudo <- function(equations,local) {
  if (equations$n!=local$n) stop("Equation/local sample sizes differ")
  if(identical(local$method,"kernel_gaussian")) return(affine_kernel_joint_pseudo(equations,local))
  PHI <- cbind(phi=equations$phi,B=equations$slope$B,b=equations$slope$b)
  out <- RdistNN::jackknife_rank_full(matrix(local$ranks,ncol=1L),PHI,
        local$rank_weights,local$deleted_rank_weights)
  base <- drop(out$theta0[,1]); deleted <- matrix(out$theta_loo[,,1],local$n,3L)
  Q <- local$n*matrix(base,local$n,3L,byrow=TRUE)-(local$n-1L)*deleted
  list(base=base,deleted=deleted,Q=Q,n=local$n)
}

affine_frozen_pseudo <- function(joint,target,application) {
  if (application=="cate") target-joint$Q[,1]
  else if (application=="slope") joint$Q[,2]*target+joint$Q[,3]
  else stop("Unknown application")
}

affine_inference <- function(equations,local,alpha=.05) {
  if (alpha!=.05) stop("The study fixes alpha=.05")
  joint <- tryCatch(affine_joint_pseudo(equations,local),error=identity)
  joint_error <- if(inherits(joint,"error")) conditionMessage(joint) else NA_character_
  if(inherits(joint,"error")) joint <- NULL
  result <- lapply(c("cate","slope"),function(app) {
    e <- equations[[app]]; w <- local$weights
    Gamma <- if(app=="cate") 1 else sum(w*e$B)
    if (!is.finite(Gamma) || Gamma<=0) stop("Nonpositive identifying coefficient")
    estimate <- -sum(w*e$b)/Gamma
    score <- e$B*estimate+e$b
    Q <- if(is.null(joint)) rep(NA_real_,local$n) else affine_frozen_pseudo(joint,estimate,app)
    variance <- c(jackknife=sum((Q-mean(Q))^2)/(local$n*(local$n-1)*Gamma^2),
        direct_moment=sum(w*score^2)/(local$N*Gamma^2),
        direct_squared=sum(w^2*score^2)/Gamma^2)
    status <- ifelse(!is.finite(variance) | variance<0,"unavailable",ifelse(variance==0,"zero","ok"))
    errors <- setNames(rep(NA_character_,length(variance)),names(variance))
    errors[status=="unavailable"] <- "Nonfinite or negative variance"
    if(!is.na(joint_error)) errors["jackknife"] <- joint_error
    SE <- sqrt(ifelse(status=="unavailable",NA_real_,variance))
    intervals <- cbind(lower=estimate-qnorm(1-alpha/2)*SE,upper=estimate+qnorm(1-alpha/2)*SE)
    list(estimate=estimate,Gamma_hat=Gamma,score_cov_scaled=variance*local$N*Gamma^2,
         variance=variance,SE=SE,intervals=intervals,variance_status=status,variance_error=errors,
         point_status=if(is.finite(estimate)) "ok" else "unavailable",pseudo_values=Q,score=score)
  }) |> setNames(c("cate","slope"))
  result$joint <- joint; result
}
