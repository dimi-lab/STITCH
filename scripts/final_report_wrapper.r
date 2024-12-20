renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--workflowpath"), type="character", default=NULL, 
              help="workflowpath. [default= %default]", metavar="character"),
  make_option(c("--data_type"), type="character", default=NULL, 
              help="data_type [default= %default]", metavar="character"),
  make_option(c("--authorname"), type="character", default=NULL, 
              help="authorname [default= %default]", metavar="character"),
  make_option(c("--samplesheet"), type="character", default=NULL, 
              help="samplesheet [default= %default]", metavar="character"),
  make_option(c("--sketch_flag"), type="character", default=NULL, 
              help="sketch_flag [default= %default]", metavar="character"),
  make_option(c("--feature_list"), type="character", default=NULL, 
              help="feature_list [default= %default]", metavar="character"),
  make_option(c("--cellcycle_correction_flag"), type="character", default=NULL, 
              help="cellcycle_correction_flag [default= %default]", metavar="character"),
  make_option(c("--vismethod"), type="character", default=NULL, 
              help="vismethod [default= %default]", metavar="character"),
  make_option(c("--seurat_obj"), type="character", default=NULL, 
              help="Path to seurat object file. [default= %default]", metavar="character"),
  make_option(c("--marker_gene_filtered"), type="character", default=NULL, 
              help="Path to marker gene filtered table. [default= %default]", metavar="character"),
  make_option(c("--diff_gene_filtered"), type="character", default=NULL, 
              help="Path to diff gene filtered table. [default= %default]", metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(opt$data_type == "scRNAseq") report_file <- file.path(opt$workflowpath,"docs/Final_report_scRNAseq.Rmd") else if(opt$data_type == "Visium") report_file <- file.path(opt$workflowpath,"docs/Final_report_Visium.Rmd") else if(opt$data_type == "VisiumHD") report_file <- file.path(opt$workflowpath,"docs/Final_report_VisiumHD.Rmd")


file.copy(report_file, './Final_report.Rmd',overwrite = T)
system(paste0('sed -i ', 's/authername/', opt$authorname, '/g ', './Final_report.Rmd'), wait = TRUE)

if(file.exists("./NO_FILE_DIFF")) opt$diff_gene_filtered <- "NA"

opt_report <- list(seurat_obj = opt$seurat_obj,
                   sketch_flag = opt$sketch_flag,
                   samplesheet = opt$samplesheet,
                   marker_gene_filtered = opt$marker_gene_filtered,
                   diff_gene_filtered = opt$diff_gene_filtered,
                   vismethod = opt$vismethod,
                   cellcycle_correction_flag = opt$cellcycle_correction_flag,
                   feature_list = opt$feature_list
                   )

params <- list(opt = opt_report)
rmarkdown::render('./Final_report.Rmd', params = params)