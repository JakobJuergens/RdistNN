#!/usr/bin/env Rscript
# Render audited production summaries. No draw generation or estimation.
argv <- commandArgs(TRUE)
if(length(argv)!=3L) stop("Usage: render_production.R ANALYSIS_DATA AUDITED_EXPORT ARTIFACT_DIR")
library(ggplot2)
library(patchwork)
data_dir <- normalizePath(argv[1]); accepted <- normalizePath(argv[2]); out <- argv[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
manifest <- jsonlite::read_json(file.path(data_dir,"analysis_manifest.json"),simplifyVector=TRUE)
stopifnot(manifest$status=="analysis_of_audited_saved_production",!manifest$estimation_performed,
  manifest$source_audit_sha256==digest::digest(file=file.path(accepted,"audit.json"),algo="sha256"))
for(i in seq_len(nrow(manifest$data_hashes))) {
  f <- manifest$data_hashes[i,]
  stopifnot(digest::digest(file=file.path(data_dir,f$file),algo="sha256")==f$sha256)
}
read <- function(root,name) read.csv(file.path(root,paste0(name,".csv")),stringsAsFactors=FALSE)
R <- read(data_dir,"raw_metrics"); D <- read(data_dir,"diagnostic_moments")
Q <- read(data_dir,"qq"); C <- read(data_dir,"paired_contrasts"); G <- read(data_dir,"gram_conditional")
P <- read(accepted,"performance"); A <- read(accepted,"paired")
V <- read(accepted,"covariance_ratios"); S <- read(accepted,"spatial_decomposition")
main <- c("full","pooled_five_fold","oracle")
key <- function(x,fields) do.call(paste,c(x[fields],sep="|"))
fields <- c("profile","cell_id","application","training_method","localization_method","variance_method")
stopifnot(!anyDuplicated(key(P,fields)),identical(key(P,fields),key(R,fields)))
P$canonical <- R$canonical
canonical <- function(x) {
  overlap <- unique(R$cell_id[R$profile=="dnn"])
  x[!(x$profile=="kernel_comparison" & x$localization_method=="dnn" & x$cell_id%in%overlap),]
}
P <- P[P$canonical & P$training_method%in%main,]
R <- R[R$canonical & R$training_method%in%main,]
A <- canonical(A); V <- canonical(V); S <- canonical(S)
S <- S[S$training_method%in%main,]; V <- V[V$training_method%in%main,]
method_names <- c(full="Full fit",pooled_five_fold="Pooled five-fold",oracle="True nuisance")
colors <- c("Full fit"="#205B80","Pooled five-fold"="#B76815","True nuisance"="#696969","Population projection"="#8B3F73",
  "DNN + jackknife"="#205B80","Kernel + jackknife"="#B76815","Kernel + conventional"="#8B3F73",
  "Jackknife"="#205B80","Direct moment"="#B76815","Squared weight"="#8B3F73")
pretty_app <- function(x) factor(ifelse(x=="cate","CATE","Residualized slope"),c("CATE","Residualized slope"))
decorate <- function(x) {
  x$app <- pretty_app(x$application)
  if("training_method"%in%names(x)) x$method <- factor(method_names[x$training_method],method_names)
  x
}
P <- decorate(P); R <- decorate(R); A$app <- pretty_app(A$application)
V <- decorate(V); S <- decorate(S); C <- decorate(C); D <- decorate(D); Q <- decorate(Q); G <- decorate(G)
theme_set(theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
  legend.position="bottom",plot.title=element_text(face="bold",size=14),
  plot.subtitle=element_text(size=10),plot.caption=element_text(size=9,hjust=0),
  strip.text=element_text(face="bold"),axis.title=element_text(size=11)))
col <- scale_color_manual(values=colors,drop=TRUE)
xlog <- function() scale_x_log10(labels=function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE))
long <- function(x,fields,labels=fields) {
  do.call(rbind,lapply(seq_along(fields),function(i) {z<-x;z$value<-x[[fields[i]]];z$measure<-factor(labels[i],levels=labels);z}))
}
plots <- receipts <- list()
save <- function(id,p,caption,width=10,height=7.5) {
  stopifnot(!id%in%names(plots))
  wrapped <- paste(strwrap(caption,width=round(width*15)),collapse="\n")
  p <- p+plot_annotation(title=gsub("_"," ",sub("^[0-9]+_","",id)),caption=wrapped)
  ggsave(file.path(out,paste0(id,".pdf")),p,device=cairo_pdf,width=width,height=height,limitsize=FALSE)
  ggsave(file.path(out,paste0(id,".png")),p,device="png",width=width,height=height,dpi=190,limitsize=FALSE)
  plots[[id]] <<- p
  receipts[[id]] <<- list(id=id,caption=caption,width_inches=width,height_inches=height)
}
base <- function(x) x[x$localization_method=="dnn" & x$order_mult==1,]
default <- base(P)
geometry <- base(R[R$training_method=="full" & R$variance_method=="jackknife" & R$application=="cate",])
g <- long(geometry,c("mean_N","mean_radius","K","mean_largest_weight"),
  c("Effective size N","Mean weighted radius","Dictionary dimension K","Largest weight"))
save("01_Local_information_and_dictionary",
  ggplot(g,aes(n,value,color=design_id))+geom_line()+geom_point()+xlog()+scale_y_log10()+
  facet_wrap(~measure,scales="free_y",ncol=2)+labs(x="Nominal sample size",y=NULL,color="Design"),
  "Full fit, default order. DNN effective size is a rank-weight quantity; physical radius also depends on reporting density. M changes at n=4,096 and 531,441.")
for(app in c("cate","slope")) {
  z <- default[default$application==app & default$design_id%in%c("C0","C1","C2") & default$variance_method=="jackknife",]
  z <- long(z,c("bias","rmse","coverage"),c("Point-target bias","RMSE","Point-target coverage"))
  save(if(app=="cate") "02_Core_CATE_performance" else "03_Core_slope_performance",
    ggplot(z,aes(n,value,color=method))+geom_line()+geom_point()+xlog()+col+
    facet_grid(measure~design_id,scales="free_y")+
    geom_hline(data=data.frame(measure=factor("Point-target coverage",levels=levels(z$measure))),aes(yintercept=.95),color="grey60")+
    labs(x="Nominal sample size",y=NULL,color=NULL),
    "10,000 common draw families per cell. Same final DNN weights for all nuisance training rules. The separate C0 endpoint changes n, M and K together.",height=8)
}
z <- P[P$design_id=="C2" & P$localization_method=="dnn" & P$variance_method=="jackknife" & P$training_method%in%c("full","oracle"),]
z$order <- factor(z$order_mult,labels=c("Order multiplier 0.4","Order multiplier 1"))
z <- long(z,c("bias_over_SD","coverage","mean_length"),c("Bias / Monte Carlo SD","Point-target coverage","Mean interval length"))
save("04_Localization_order_tradeoff",ggplot(z,aes(n,value,color=method,linetype=order))+
  geom_line()+geom_point()+xlog()+col+facet_grid(measure~app,scales="free_y")+labs(x="Nominal sample size",y=NULL,color=NULL,linetype=NULL),
  "C2 has reporting curvature represented by q5. Higher order concentrates weights: lower observed spatial drift here accompanies wider intervals. This is not a monotonicity theorem.",height=8)
z <- base(A[A$design_id%in%c("C0","C1","C2"),]); z$gap <- factor(z$comparison,
  levels=c("full_vs_CF","full_vs_oracle"),labels=c("Full minus five-fold","Full minus true nuisance"))
save("05_Training_and_oracle_gaps",ggplot(z,aes(n,q95_absolute_Omega_units,color=gap))+
  geom_line()+geom_point()+xlog()+scale_y_log10()+facet_grid(app~design_id)+
  labs(x="Nominal sample size",y="95th percentile of absolute normalized gap / point-limit SD",color=NULL),
  "Normalization uses sqrt(N) and the exact point-limit target covariance. Close full/five-fold estimates can share a persistent true-nuisance gap; these are paired gaps, not separate SDs.")
z <- C[C$comparison=="full_minus_pooled_five_fold" & C$localization_method=="dnn" & C$order_mult==1,]
z <- long(z,c("coverage_difference","length_difference","MSE_difference"),c("Coverage difference","Interval-length difference","MSE difference"))
z$uncertainty <- ifelse(z$measure=="Coverage difference",z$coverage_difference_mcse,
  ifelse(z$measure=="Interval-length difference",z$length_difference_mcse,z$MSE_difference_mcse))
save("06_Paired_full_minus_five_fold",ggplot(z,aes(n,value,color=design_id))+
  geom_hline(yintercept=0,color="grey60")+geom_errorbar(aes(ymin=value-1.96*uncertainty,ymax=value+1.96*uncertainty),width=0)+
  geom_line()+geom_point()+xlog()+facet_grid(measure~app,scales="free_y")+labs(x="Nominal sample size",y=NULL,color="Design"),
  "Left minus right, 10,000 aligned pairs; bars are +/- 1.96 paired Monte Carlo SE, not simultaneous bands. Zero observed coverage difference does not imply identical interval decisions.",height=8)
z <- base(V[V$design_id=="C0" & V$training_method=="full",])
z$panel <- paste(z$scale,z$reference,sep=" / ")
save("07_Covariance_reference_comparisons",ggplot(z,aes(n,ratio_of_means,color=panel))+
  geom_hline(yintercept=1,color="grey60")+geom_line()+geom_point()+xlog()+facet_wrap(~app)+
  labs(x="Nominal sample size",y="Ratio of mean covariances",color="Scale / reference"),
  "Point-oracle array, oracle fitted-root array and asymptotic point constant are different references. Score-to-target conversion uses the fitted or oracle identifying coefficient as appropriate.")
z <- default[default$design_id=="C0",]
z$route <- factor(z$variance_method,levels=c("jackknife","direct_moment","direct_squared"),
  labels=c("Jackknife","Direct moment","Squared weight"))
z <- long(z,c("rms_SE_over_SD","coverage"),c("RMS standard error / Monte Carlo SD","Point-target coverage"))
save("08_DNN_variance_routes",ggplot(z,aes(n,value,color=route,linetype=method))+
  geom_line()+geom_point()+xlog()+col+facet_grid(measure~app,scales="free_y")+
  labs(x="Nominal sample size",y=NULL,color=NULL,linetype=NULL),
  "Variance routes share each point estimate. RMS SE/SD measures second-moment calibration; point-target coverage also depends on centering and studentized shape. No universal jackknife advantage is asserted.")
dd <- base(D[D$design_id=="C0" & D$training_method=="full" & D$observed>0,])
z <- dd[dd$field%in%paste0("covariance_",c("transfer","component","target","cross"),"_ratio_point_oracle"),]
z$component <- sub("_ratio_point_oracle","",sub("covariance_","",z$field))
save("09_Signed_covariance_transfer_terms",ggplot(z,aes(n,mean,color=component))+
  geom_hline(yintercept=0,color="grey60")+geom_line()+geom_point()+xlog()+facet_wrap(~app)+
  labs(x="Nominal sample size",y="Mean term / draw-specific point-oracle score covariance",color="Transfer term"),
  "Transfer square T=component+target+cross. Cross is signed and already includes the factor two. T is not fitted-minus-oracle covariance. Target movement cancels for CATE, remains for slope.")
z <- dd[dd$application=="slope" & dd$field%in%paste0("covariance_",c("target","cross"),"_ratio_point_oracle"),]
save("10_Slope_target_movement_detail",ggplot(z,aes(n,mean,color=field))+
  geom_hline(yintercept=0,color="grey60")+geom_line()+geom_point()+xlog()+
  scale_color_discrete(labels=c("Signed cross term","Target movement square"))+
  labs(x="Nominal sample size",y="Mean term / point-oracle score covariance",color=NULL),
  "Movement is small in these cells but nonzero and nonmonotone. Its small numerical magnitude does not justify omitting it from the general varying-coefficient proof.",height=5)
z <- dd[dd$field%in%paste0("covariance_",c("transfer","oracle_cross"),"_ratio_point_oracle"),]
save("11_Transfer_and_covariance_inflation",ggplot(z,aes(n,mean,color=field))+
  geom_hline(yintercept=0,color="grey60")+geom_line()+geom_point()+xlog()+facet_wrap(~app)+
  scale_color_discrete(labels=c("Oracle cross: F - O - T","Transfer square T"))+
  labs(x="Nominal sample size",y="Mean ratio to point-oracle score covariance",color=NULL),
  "F/O = 1 + T/O + oracle-cross/O draw by draw. Oracle cross is obtained by subtraction of saved components; it is not independently reconstructed as an inner product.")
case_label <- function(z) ifelse(z$design_id=="C0" & z$n==531441,"C0, M3 endpoint",
  paste0(z$design_id,", n=",format(z$n,big.mark=",",scientific=FALSE,trim=TRUE),", order ",z$order_mult))
selected <- function(z) (z$design_id=="C0" & z$n%in%c(512,65536,531441)) |
  (z$design_id=="C1" & z$n==2048) | (z$design_id=="C2" & z$n==512 & z$order_mult==.4) |
  (z$design_id=="C2S" & z$n==2048)
for(app in c("cate","slope")) for(def in c("point_t","centered_t")) {
  z <- Q[selected(Q) & Q$application==app & Q$variance_method=="jackknife" & Q$localization_method=="dnn" & Q$definition==def,]
  z$case <- case_label(z)
  id <- if(app=="cate") if(def=="point_t") "12_CATE_point_student_QQ" else "13_CATE_empirical_center_QQ" else
    if(def=="point_t") "14_Slope_point_student_QQ" else "15_Slope_empirical_center_QQ"
  save(id,ggplot(z,aes(normal_quantile,empirical_quantile,color=method))+
    geom_abline(slope=1,intercept=0,color="grey60")+geom_line()+col+facet_wrap(~case,ncol=3)+
    labs(x="Standard normal quantile",y=if(def=="point_t") "(Estimate - point target) / fitted SE" else "(Estimate - Monte Carlo mean) / fitted SE",color=NULL),
    paste0("Genuine type-8 quantiles at 101 probabilities, 10,000 observations per curve; no quantile bands. ",
      if(def=="centered_t") "Monte Carlo centering is an unavailable diagnostic, not a proposed bias correction." else "All methods use their own frozen-score SE; point centering retains both nuisance and localization shifts."))
}
z <- default[default$variance_method=="jackknife" & default$training_method%in%c("full","oracle") & default$design_id%in%c("C0","C1","C2"),]
z <- long(z,c("student_lower_tail","student_upper_tail"),c("Lower student tail","Upper student tail"))
save("16_Point_student_tail_balance",ggplot(z,aes(n,value,color=method,linetype=measure))+
  geom_hline(yintercept=.025,color="grey60")+geom_line()+geom_point()+xlog()+col+
  facet_grid(app~design_id)+labs(x="Nominal sample size",y="Fraction beyond +/- 1.96",color=NULL,linetype=NULL),
  "Point student tails, 10,000 observations. Balanced excess tails differ from one-sided centering failure. Normal reference is .025 per tail.")
kernel <- P[P$profile%in%c("kernel_comparison","settings_comparison") & P$order_mult==1 & P$training_method=="full" &
  ((P$localization_method=="dnn" & P$variance_method=="jackknife") |
    (P$localization_method=="kernel_gaussian" & P$variance_method%in%c("jackknife","direct_squared"))),]
# Reintroduce only the canonical DNN endpoints present in the kernel comparison.
cells <- unique(kernel$cell_id)
add <- P[P$profile=="dnn" & P$cell_id%in%cells & P$training_method=="full" & P$variance_method=="jackknife",]
kernel <- rbind(kernel,add)
kernel$procedure <- ifelse(kernel$localization_method=="dnn","DNN + jackknife",
  ifelse(kernel$variance_method=="jackknife","Kernel + jackknife","Kernel + conventional"))
stopifnot(!anyDuplicated(key(kernel,c("cell_id","application","procedure"))))
for(app in c("cate","slope")) {
  z <- long(kernel[kernel$application==app,],c("coverage","mean_length","rmse"),c("Point-target coverage","Mean interval length","RMSE"))
  save(if(app=="cate") "17_Feasible_kernel_comparison_CATE" else "18_Feasible_kernel_comparison_slope",
    ggplot(z,aes(n,value,color=procedure))+geom_line()+geom_point()+xlog()+col+
    facet_grid(measure~design_id,scales="free_y")+labs(x="Nominal sample size",y=NULL,color=NULL),
    "Kernel bandwidth matches population weight concentration to DNN, not realized ESS or bias. Same nuisance predictions. Kernel JK/conventional changes only the variance route; DNN/kernel also changes localization.",width=12,height=8)
}
z <- C[C$training_method=="full" & C$order_mult==1 &
  (C$comparison=="dnn_jk_minus_kernel_sandwich" |
   (C$comparison=="jackknife_minus_direct_squared" & C$localization_method=="kernel_gaussian")),]
z$comparison <- ifelse(z$comparison=="dnn_jk_minus_kernel_sandwich","DNN JK - kernel conventional","Kernel JK - kernel conventional")
z <- long(z,c("coverage_difference","length_difference"),c("Coverage difference","Interval-length difference"))
z$uncertainty <- ifelse(z$measure=="Coverage difference",z$coverage_difference_mcse,z$length_difference_mcse)
save("19_Paired_kernel_comparison_costs",ggplot(z,aes(n,value,color=design_id,shape=comparison))+
  geom_hline(yintercept=0,color="grey60")+geom_errorbar(aes(ymin=value-1.96*uncertainty,ymax=value+1.96*uncertainty),width=0)+
  geom_point(position=position_dodge(width=.02))+xlog()+facet_grid(measure~app,scales="free_y")+
  labs(x="Nominal sample size",y=NULL,color="Design",shape=NULL),
  "Left minus right, default order, 10,000 pairs. Bars are +/- 1.96 paired MCSE. Coverage gains have interval-length costs; cross-procedure contrasts do not isolate variance estimation.",height=8)
z <- C[C$comparison%in%c("C0H_minus_baseline","C2S_minus_baseline") & C$training_method=="full" &
  C$localization_method=="dnn" & C$left_route=="jackknife",]
z <- long(z,c("coverage_difference","MSE_difference"),c("Coverage difference","MSE difference"))
z$uncertainty <- ifelse(z$measure=="Coverage difference",z$coverage_difference_mcse,z$MSE_difference_mcse)
save("20_Paired_settings_changes",ggplot(z,aes(n,value,color=design_id))+
  geom_hline(yintercept=0,color="grey60")+geom_errorbar(aes(ymin=value-1.96*uncertainty,ymax=value+1.96*uncertainty),width=0)+
  geom_line()+geom_point()+xlog()+facet_grid(measure~app,scales="free_y")+labs(x="Nominal sample size",y=NULL,color="Setting minus baseline"),
  "C0H - C0 and C2S - C2, common latent draws. Sparse density changes physical radius at the same DNN rank-weight ESS. Bars show paired Monte Carlo uncertainty.")
z <- D[D$design_id%in%c("C0","C1") & D$localization_method=="dnn" & D$order_mult==1 & D$observed>0 &
  ((D$application=="cate" & D$field=="point_cate_prediction_MSE_sum") |
    (D$application=="slope" & D$field=="point_slope_prediction_MSE_sum")) &
  D$training_method%in%c("full","pooled_five_fold","population_projection"),]
z$method <- ifelse(z$training_method=="population_projection","Population projection",as.character(z$method))
save("21_Fresh_point_prediction_risk",ggplot(z,aes(n,mean,color=method))+
  geom_line()+geom_point()+xlog()+scale_y_log10()+col+facet_grid(app~design_id)+
  labs(x="Nominal sample size",y="Mean sum of relevant fresh point prediction MSEs",color=NULL),
  "CATE sums propensity and two outcome-mean risks; slope sums propensity and marginal-outcome risks. Fitted risks observe 10,000 draws; projection controls use 1,000, or 100 at M3. These are point risks, not contract suprema.")
z <- S[S$training_method=="full" & S$localization_method=="dnn" & S$design_id%in%c("C0","C1","C2","C2S"),]
z$case <- paste0(z$design_id,", order ",z$order_mult)
z <- long(z,c("mean_drift","mean_remainder"),c("Realized weighted-truth drift","Estimate minus realized weighted truth"))
save("22_Point_centering_decomposition",ggplot(z,aes(n,value,color=measure))+
  geom_line()+geom_point()+xlog()+facet_grid(app~case,scales="free_y")+
  labs(x="Nominal sample size",y="Mean target-scale component",color=NULL),
  "Exact sample decomposition. Remainder includes oracle outcome noise and nuisance effects; it is not pure orthogonal mean bias. Slope realized weighted truth uses residual-variance weights.",width=12)
z <- G[G$design_id%in%c("C0","C1") & G$n==4096,]
z$cohort <- ifelse(z$gram_event,"Gram event holds (9,979)","Gram event fails (21)")
save("23_Gram_event_cohort_coverage",ggplot(z,aes(cohort,coverage,color=design_id))+
  geom_hline(yintercept=.95,color="grey60")+
  geom_errorbar(aes(ymin=coverage_lower,ymax=coverage_upper),position=position_dodge(width=.3),width=.12)+
  geom_point(position=position_dodge(width=.3))+facet_wrap(~app)+
  labs(x=NULL,y="Point-target coverage",color="Design"),
  "C0 and C1 share the same 21 failing first-stage Gram draws; observations remain in all production summaries. Exact binomial intervals are descriptive within these selected cohorts, not a causal event comparison.")
# One page per figure; standalone PDF/PNG pairs remain the reusable artifacts.
cairo_pdf(file.path(out,"production_atlas.pdf"),width=12,height=8.5,onefile=TRUE)
for(p in plots) print(p)
dev.off()
files <- list.files(out,pattern="\\.(pdf|png)$",full.names=TRUE)
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE))
receipt <- list(status="rendered_saved_production_analysis",simulation_performed=FALSE,
  renderer_sha256=digest::digest(file=script,algo="sha256"),
  analysis_manifest_sha256=digest::digest(file=file.path(data_dir,"analysis_manifest.json"),algo="sha256"),
  accepted_audit_sha256=manifest$source_audit_sha256,figures=receipts,
  output_hashes=lapply(files,function(f) list(file=basename(f),sha256=digest::digest(file=f,algo="sha256"))))
jsonlite::write_json(receipt,file.path(out,"manifest.json"),pretty=TRUE,auto_unbox=TRUE,digits=17)
index <- c("# Affine production simulation atlas","", "Generated from saved, audited 10,000-repetition production. No estimation was rerun.","",
  "The atlas is an analysis artifact for manuscript selection. It is not a replacement manuscript or a claim of theorem verification.","",
  "[Combined PDF atlas](production_atlas.pdf)","", "| Figure | Interpretation |","| --- | --- |")
for(id in names(receipts)) index <- c(index,paste0("| [",id,"](",id,".pdf) / [PNG](",id,".png) | ",receipts[[id]]$caption," |"))
writeLines(index,file.path(out,"README.md"))
cat(length(plots),"figure pairs and",length(plots),"atlas pages generated from saved production.\n")
