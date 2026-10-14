#!/usr/bin/env Rscript
# Publish validated portable SIPP tables/PDFs; no nuisance fitting or RDS reads.
.script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1L])
.dir <- dirname(normalizePath(.script))
repo <- normalizePath(file.path(.dir,"../../.."))
source(file.path(.dir,"R","manuscript_export.R"))
args <- commandArgs(TRUE)
options <- list(input_dir=file.path(repo,"markdown_notes/affine_scores/sipp_development/2026-10-12_code_overhaul"),
  historical_dir=file.path(repo,"markdown_notes/affine_scores/sipp_development/2026-10-12_overhaul/evidence"),
  output_dir=file.path(repo,"Artifacts/affine_scores/sipp_application_20261013"),
  table_dir=file.path(repo,"Tables/Applications/SIPP_401k/affine_overhaul_20261013"),check_only=FALSE)
for(arg in args) {
  if(arg=="--check"){options$check_only<-TRUE;next}
  key<-sub("=.*$","",arg); value<-sub("^[^=]*=","",arg)
  field<-switch(key,"--input-dir"="input_dir","--historical-dir"="historical_dir",
    "--output-dir"="output_dir","--table-dir"="table_dir",stop("Unknown argument: ",arg,call.=FALSE))
  options[[field]]<-value
}
do.call(render_sipp_401k_manuscript,c(options,list(repo_root=repo)))
cat(if(options$check_only)"Checked" else "Published","SIPP manuscript tables and three verified PDF displays.\n")
