# One model authority for estimators, controls, and fresh-law diagnostics.
affine_reporting_quantile <- function(u,design_id="C0") {
  if(any(!is.finite(u)) || any(u<0 | u>1)) stop("Invalid reporting quantiles")
  if(!affine_model_profile(design_id)$sparse) return(u)
  ifelse(u<3/8,u/1.5,ifelse(u<=5/8,.25+2*(u-3/8),.75+(u-5/8)/1.5))
}

affine_reporting_density <- function(u,design_id="C0") {
  sparse <- affine_model_profile(design_id)$sparse
  ifelse(u<0 | u>1,0,if(sparse) .5+as.numeric(u<.25 | u>.75) else 1)
}

affine_reporting_gram <- function(q_dim,design_id="C0") {
  if(!q_dim%in%c(3,5)) stop("Invalid fixed reporting basis")
  diag(c(1,rep(if(affine_model_profile(design_id)$sparse) 11/8 else 1,2),rep(1,q_dim-3)))
}

affine_conditional_noise_variance <- function(x1,design_id="C0") {
  if(affine_model_profile(design_id)$heteroskedastic) (3*x1^2+9*x1+7)/784 else rep(1/64,length(x1))
}

affine_noise_cell_variance <- function(u,bins,M,design_id="C0") {
  if(!affine_model_profile(design_id)$heteroskedastic) return(rep(1/64,nrow(bins)))
  lo <- (bins[,1]-1)/M; hi <- bins[,1]/M
  3/784*((1+u[1])^2+(1+u[1])*(lo+hi)+(lo^2+lo*hi+hi^2)/3)
}

affine_reporting_ball <- function(r,design_id="C0") {
  if(any(!is.finite(r)) || any(r<0)) stop("Invalid reporting radius")
  rectangle <- function(a,b) {
    end <- pmin(a,r); flat <- pmin(end,sqrt(pmax(r^2-b^2,0)))
    primitive <- function(t) .5*(t*sqrt(pmax(r^2-t^2,0))+r^2*asin(pmin(t/pmax(r,.Machine$double.xmin),1)))
    4*(b*flat+primitive(end)-primitive(flat))
  }
  F <- if(affine_model_profile(design_id)$sparse)
    2.25*rectangle(.5,.5)-3*rectangle(.25,.5)+rectangle(.25,.25) else rectangle(.5,.5)
  pmin(1,pmax(0,F))
}

affine_unit_matrix <- function(x, columns, name) {
  if (!is.matrix(x) || !is.numeric(x) || ncol(x)!=columns || any(!is.finite(x)) ||
      any(x<0 | x>1)) stop(paste(name,"must be a finite unit-cube matrix"))
  invisible(TRUE)
}

affine_latent <- function(n, seed) {
  if (n < 1 || n != floor(n)) stop("Invalid latent size")
  affine_with_seed(seed, matrix(runif(n*9L),n,9L,byrow=TRUE))
}

affine_truth <- function(X, Z, design_id = "C0") {
  affine_unit_matrix(X,2,"X"); affine_unit_matrix(Z,4,"Z")
  if (nrow(X)!=nrow(Z)) stop("Invalid model inputs")
  profile <- affine_model_profile(design_id)
  S <- sin(2*pi*Z)
  propensity <- .5 + S[,4]/8
  baseline <- (S[,1]+X[,1]*S[,2]+X[,2]*S[,3])/8
  if (profile$base=="C1") baseline <- baseline + S[,4]/8
  theta <- rowSums(X)/8
  if (profile$base=="C2") theta <- theta + (2-cos(8*pi*(X[,1]-.5))-cos(8*pi*(X[,2]-.5)))/16
  list(pi=propensity,mu0=baseline,mu1=baseline+theta,m=baseline+propensity*theta,
       theta=theta,noise_scale=if(profile$heteroskedastic) sqrt(3)/28*(1+X[,1]+Z[,1]) else rep(1/8,nrow(X)))
}

affine_data <- function(latent, design_id = "C0") {
  affine_unit_matrix(latent,9,"latent")
  X <- affine_reporting_quantile(latent[,1:2,drop=FALSE],design_id); Z <- latent[,3:6,drop=FALSE]
  truth <- affine_truth(X,Z,design_id)
  T <- as.integer(latent[,7] < truth$pi)
  eps0 <- ifelse(latent[,8]<.5,-truth$noise_scale,truth$noise_scale)
  eps1 <- ifelse(latent[,9]<.5,-truth$noise_scale,truth$noise_scale)
  Y <- truth$mu0+T*truth$theta+(1-T)*eps0+T*eps1
  list(X=X,Z=Z,T=T,Y=Y,truth=truth,n=nrow(X),design_id=design_id)
}

affine_basis <- function(X, q_dim = 3L) {
  affine_unit_matrix(X,2,"X")
  if (!q_dim %in% c(3,5)) stop("Only fixed q3/q5 factors are supported")
  Q <- cbind(1,sqrt(12)*(X[,1]-.5),sqrt(12)*(X[,2]-.5))
  if (q_dim==5) Q <- cbind(Q,sqrt(2)*cos(8*pi*X[,1]),sqrt(2)*cos(8*pi*X[,2]))
  Q
}

affine_cells <- function(Z, M) {
  affine_unit_matrix(Z,4,"Z")
  if (M<1 || M!=floor(M)) stop("Invalid histogram resolution")
  bins <- pmin(floor(M*Z),M-1)
  as.integer(1+drop(bins %*% M^(0:3)))
}

affine_sine_cells <- function(M) {
  lo <- (0:(M-1))/M; hi <- (1:M)/M
  mean <- M*(cos(2*pi*lo)-cos(2*pi*hi))/(2*pi)
  second <- .5-M*(sin(4*pi*hi)-sin(4*pi*lo))/(8*pi)
  data.frame(mean=mean,second=second,variance=pmax(0,second-mean^2))
}

affine_projection_bias <- function(M) {
  S <- affine_sine_cells(M); A <- S$variance/64; p <- .5+S$mean/8
  list(A=mean(A),cate=mean(A/(p*(1-p))),slope=mean(A)/(31/128+mean(A)))
}

affine_cell_moments <- function(M) {
  bins <- as.matrix(expand.grid(rep(list(seq_len(M)),4)))
  S <- affine_sine_cells(M)
  list(bins=bins,sine_mean=matrix(S$mean[bins],nrow(bins),4),
       sine_variance=matrix(S$variance[bins],nrow(bins),4))
}

affine_population_coefficients <- function(schedule, design_id = schedule$design_id, ridge = FALSE) {
  profile <- affine_model_profile(design_id)
  C <- affine_cell_moments(schedule$M); S <- C$sine_mean; J <- schedule$J; qd <- schedule$q_dim
  p <- .5+S[,4]/8
  c0 <- cbind((S[,1]+(S[,2]+S[,3])/2)/8,S[,2]/(8*sqrt(12)),S[,3]/(8*sqrt(12)))
  if (profile$base=="C1") c0[,1] <- c0[,1]+S[,4]/8
  if (qd==5) c0 <- cbind(c0,0,0)
  d <- c(1/8,1/(8*sqrt(12)),1/(8*sqrt(12)))
  if (profile$base=="C2") {
    d[1] <- 1/4
    if (qd==5) d <- c(d,-1/(16*sqrt(2)),-1/(16*sqrt(2)))
  } else if (qd==5) d <- c(d,0,0)
  dm <- matrix(d,J,qd,byrow=TRUE)
  out <- list(pi=cbind(p,matrix(0,J,qd-1)),mu0=c0,mu1=c0+dm,m=c0+p*dm)
  if (profile$base=="C1") {
    A <- C$sine_variance[,4]/64
    out$mu0[,1] <- out$mu0[,1]-A/(1-p)
    out$mu1[,1] <- out$mu1[,1]+A/p
  }
  if (ridge) {
    lambda <- schedule$lambda
    if(profile$sparse) {
      G <- diag(affine_reporting_gram(qd,design_id))
      for(j in names(out)) {
        arm <- switch(j,mu0=1-p,mu1=p,rep(1,J))
        for(k in seq_len(qd)) out[[j]][,k] <- out[[j]][,k]*arm*G[k]/(arm*G[k]+lambda)
      }
    } else {
      out$pi <- out$pi/(1+lambda); out$m <- out$m/(1+lambda)
      out$mu0 <- out$mu0*((1-p)/(1-p+lambda)); out$mu1 <- out$mu1*(p/(p+lambda))
    }
  }
  out
}
