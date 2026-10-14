#!/usr/bin/env Rscript
# Focused independent checks of saved-data analysis and rendered receipts.
args <- commandArgs(TRUE)
if(length(args)!=3L) stop("Usage: verify_production_analysis.R ANALYSIS_DATA ACCEPTED_EXPORT ARTIFACT_DIR")
d <- args[1]; a <- args[2]; out <- args[3]
read <- function(root,f) read.csv(file.path(root,paste0(f,".csv")),stringsAsFactors=FALSE)
M <- jsonlite::read_json(file.path(d,"analysis_manifest.json"),simplifyVector=TRUE)
stopifnot(!M$estimation_performed,M$source_audit_sha256==digest::digest(file=file.path(a,"audit.json"),algo="sha256"))
for(i in seq_len(nrow(M$data_hashes))) {
  h <- M$data_hashes[i,];stopifnot(digest::digest(file=file.path(d,h$file),algo="sha256")==h$sha256)
}
R <- read(d,"raw_metrics"); P <- read(a,"performance"); Q <- read(d,"qq")
C <- read(d,"paired_contrasts"); D <- read(d,"diagnostic_moments"); G <- read(d,"gram_conditional")
key <- function(x,f) do.call(paste,c(x[f],sep="|"))
keys <- c("profile","cell_id","application","training_method","localization_method","variance_method")
stopifnot(nrow(R)==918,nrow(Q)==179982,nrow(C)==816,nrow(G)==44,
  identical(key(R,keys),key(P,keys)),!anyDuplicated(key(R,keys)))
for(pair in list(c("bias","bias"),c("sd","sd"),c("coverage_point","coverage"),
 c("mean_N","mean_N"),c("mean_radius","mean_weighted_radius"),c("mean_largest_weight","mean_largest_weight"),
 c("point_t_mean","student_mean"),c("point_t_sd","student_sd")))
  stopifnot(isTRUE(all.equal(R[[pair[1]]],P[[pair[2]]],tolerance=1e-13)))
stopifnot(all(C$paired==10000),all(C$point_pairs==10000),all(C$excluded==0),all(C$excluded_points==0))
stopifnot(max(abs(C$coverage_difference-(C$left_coverage-C$right_coverage)))<1e-14,
  max(abs(C$coverage_difference-(C$coverage_gains-C$coverage_losses)/C$paired))<1e-14,
  all(C$coverage_discordances==C$coverage_gains+C$coverage_losses))
v <- (C$coverage_discordances-C$paired*C$coverage_difference^2)/(C$paired-1)
stopifnot(max(abs(sqrt(v/C$paired)-C$coverage_difference_mcse))<1e-13)
qkeys <- c(keys,"definition")
curve <- split(seq_len(nrow(Q)),key(Q,qkeys))
stopifnot(length(curve)==1782,all(Q$observed==10000),all(vapply(curve,function(ix)
  length(ix)==101 && all(diff(Q$probability[ix])>0) && all(diff(Q$empirical_quantile[ix])>=0),logical(1))))
t <- Q[Q$definition=="point_t" & Q$probability%in%c(.025,.975),]
ix <- match(key(t,keys),key(P,keys));stopifnot(!anyNA(ix))
expected <- ifelse(t$probability==.025,P$student_q025[ix],P$student_q975[ix])
stopifnot(max(abs(t$empirical_quantile-expected))<1e-13)
f <- c("profile","cell_id","application","training_method","localization_method")
covfields <- paste0("covariance_",c("transfer","component","target","cross"))
wide <- function(field) {z<-D[D$field==field & D$observed>0,];setNames(z$mean,key(z,f))}
v <- lapply(covfields,wide);common<-Reduce(intersect,lapply(v,names))
stopifnot(max(abs(v[[1]][common]-v[[2]][common]-v[[3]][common]-v[[4]][common]))<1e-13)
stopifnot(max(abs(D$maximum[D$application=="cate" & D$field=="covariance_target" & D$observed>0]))<1e-25)
for(ix in split(seq_len(nrow(G)),key(G,f))) {
  z <- G[ix,];j<-match(key(z[1,],f),key(R[R$variance_method=="jackknife",],f))
  r<-R[R$variance_method=="jackknife",][j,]
  stopifnot(sum(z$observed)==10000,abs(weighted.mean(z$coverage,z$observed)-r$coverage_point)<1e-13,
    abs(weighted.mean(z$bias,z$observed)-r$bias)<1e-13)
}
proj <- read(d,"projection_references"); z<-proj[proj$design_id=="C0",]
stopifnot(max(abs(z$cate_MSE_sum-z$sine_residual_variance/16))<1e-13,
  z$cate_MSE_sum[z$M==3]>z$cate_MSE_sum[z$M==2])
F <- jsonlite::read_json(file.path(out,"manifest.json"),simplifyVector=TRUE)
stopifnot(!F$simulation_performed,length(F$figures)==23,nrow(F$output_hashes)==47,
  F$analysis_manifest_sha256==digest::digest(file=file.path(d,"analysis_manifest.json"),algo="sha256"))
for(i in seq_len(nrow(F$output_hashes))) {
  h<-F$output_hashes[i,];stopifnot(digest::digest(file=file.path(out,h$file),algo="sha256")==h$sha256)
}
cat("PASS: six table hashes, 918 shared rows, 1,782 QQ curves, 816 paired contrasts, covariance identity, 44 Gram cohorts, analytic projection and 47 rendered-file hashes.\n")
