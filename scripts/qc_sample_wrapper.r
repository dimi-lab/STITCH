renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--workflowpath"), type="character", default=NULL, 
              help="workflow path. [default= %default]", metavar="character"),
  make_option(c("--data_type"), type="character", default=NULL, 
              help="data type. [default= %default]", metavar="character"),
  make_option(c("--binsize"), type="integer", default=NULL, 
              help="bin size to read in. [default= %default]", metavar="integer"),
  make_option(c("--ambient_RNA_removal_flag"), type="character", default=NULL, 
              help="ambient RNA removal flag. [default= %default]", metavar="character"),
  make_option(c("--doublet_removal_flag"), type="character", default=NULL, 
              help="doublet removal flag. [default= %default]", metavar="character"),
  make_option(c("--adaptive_cutoff_flag"), type="character", default=NULL, 
              help="adaptive cutoff flag. [default= %default]", metavar="character"),
  make_option(c("--mt_cutoff"), type="character", default=NULL, 
              help="mt cutoff. [default= %default]", metavar="character"),
  make_option(c("--hb_cutoff"), type="character", default=NULL, 
              help="hb cutoff. [default= %default]", metavar="character"),
  make_option(c("--nFeature_cutoff"), type="character", default=NULL, 
              help="nFeature cutoff. [default= %default]", metavar="character"),
  make_option(c("--nCount_cutoff"), type="character", default=NULL, 
              help="nCount cutoff. [default= %default]", metavar="character"),
  make_option(c("--nCell_cutoff"), type="character", default=NULL, 
              help="nCell cutoff. [default= %default]", metavar="character"),
  make_option(c("--sampleid"), type="character", default=NULL, 
              help="sample id. [default= %default]", metavar="integer"),
  make_option(c("--condition"), type="character", default=NULL, 
            help="condition. [default= %default]", metavar="integer"),
  make_option(c("--secondary_output"), type="character", default=NULL, 
            help="Path to secondary output. [default= %default]", metavar="integer")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

options(stringsAsFactors = FALSE)

for(i in c("mt_cutoff", "hb_cutoff", "nFeature_cutoff", "nCount_cutoff","nCell_cutoff", "binsize")) opt[[i]] <- as.numeric(opt[[i]])

params <- list(opt = opt)
message("generate QC report")

if(opt$data_type == 'scRNAseq') {
  file.copy(file.path(opt$workflowpath, "/docs/QC_report_scRNAseq.Rmd"), "./")
  rmarkdown::render("./QC_report_scRNAseq.Rmd", params = params, output_file = "QC_report.html")
  file.remove("./QC_report_scRNAseq.Rmd")
} else if(opt$data_type == 'Visium'){
  file.copy(file.path(opt$workflowpath, "/docs/QC_report_Visium.Rmd"), "./")
  rmarkdown::render("./QC_report_Visium.Rmd", params = params, output_file = "QC_report.html")
  file.remove("./QC_report_Visium.Rmd")
} else if(opt$data_type == 'VisiumHD'){
  file.copy(file.path(opt$workflowpath, "/docs/QC_report_VisiumHD.Rmd"), "./")
  rmarkdown::render("./QC_report_VisiumHD.Rmd", params = params, output_file = "QC_report.html")
  file.remove("./QC_report_VisiumHD.Rmd")
}