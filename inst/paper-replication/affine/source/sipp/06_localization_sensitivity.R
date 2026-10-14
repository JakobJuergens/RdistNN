# Broader/narrower localization from saved SIPP predictions; no refitting.
.script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1L])
.dir<-dirname(normalizePath(.script))
source(file.path(.dir,"R","estimation_helpers.R"))
source(file.path(.dir,"R","localization_sensitivity.R"))
options<-list(fit_cache=NULL,output_dir=NULL,orders=NULL,allow_validation=FALSE)
for(arg in commandArgs(trailingOnly=TRUE)) {
  if(arg=="--allow-validation"){options$allow_validation<-TRUE;next}
  key<-sub("=.*$","",arg); value<-sub("^[^=]*=","",arg)
  field<-switch(key,"--fit-cache"="fit_cache","--output-dir"="output_dir","--orders"="orders",
    stop("Unknown localization argument: ",arg,call.=FALSE))
  options[[field]]<-if(field=="orders")suppressWarnings(as.numeric(strsplit(value,",",fixed=TRUE)[[1L]])) else value
}
if(is.null(options$fit_cache)||is.null(options$output_dir)||is.null(options$orders))
  stop("Provide --fit-cache, --orders and a new --output-dir.",call.=FALSE)
export_sipp_401k_localization(readRDS(options$fit_cache),options$orders,options$output_dir,
                            options$allow_validation)
cat("Exported saved-fit localization sensitivity:",options$output_dir,"\n")
