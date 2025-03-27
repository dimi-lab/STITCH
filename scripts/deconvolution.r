renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--reference"), type="character", default=NULL, 
              help="Path to reference seurat object in .rds format. [default= %default]", metavar="character"),
  make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample id. [default= %default]", metavar="character"),
  make_option(c("--reference_assay"), type="character", default=NULL, 
              help="Reference Assay to use to perform deconvolution. [default= %default]", metavar="character"),
  make_option(c("--query_assay"), type="character", default=NULL, 
              help="Qeury assay to use to perform deconvolution [default= %default]", metavar="character"),
  make_option(c("--refdata"), type="character", default=NULL, 
              help="Character indicating reference data, e.g. cell type, to use for deconvolution. Must exist in the metadata of the reference seurat object. [default= %default]", metavar="character"),
  make_option(c("--doublet_mode"), type="character", default=NULL, 
              help="doublet mode for RCTD, can be either doublet, multi or full. See help page for run.RCTD for details. [default= %default]", metavar="character"),
  make_option(c("--parallel_strategy"), type="character", default=NULL, 
              help="Parallel strategy for future. See help page for plan for details. [default= %default]", metavar="character"),
  make_option(c("--nworkers"), type="integer", default=NULL, 
              help="Number of workers/cpus used for future. [default= %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(Seurat)
library(spacexr)
library(future)
library(doFuture)
registerDoFuture()
plan(opt$parallel_strategy, workers = as.integer(opt$nworkers))
options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)

remove_object <- function(object_name){
  for(i in object_name) assign(i, NULL,envir = .GlobalEnv)
  rm(list = object_name, envir = .GlobalEnv)
  invisible(gc(verbose = FALSE))
}

reference <- readRDS(opt$reference)
query <- readRDS(paste0(opt$sampleid, ".rds"))

reference <- Reference(reference[[opt$reference_assay]]$counts, as.factor(reference@meta.data[[opt$refdata]]))

bulk_spatial <- SpatialRNA(GetTissueCoordinates(query)[,1:2], query[[opt$query_assay]]$counts)
myRCTD <- create.RCTD(bulk_spatial, reference, max_cores = 10, CELL_MIN_INSTANCE=min(table(reference@cell_types)), UMI_min = 0, counts_MIN = 0)
myRCTD <- run.RCTD(myRCTD, doublet_mode = opt$doublet_mode)

message("saving RCTD object")
saveRDS(myRCTD@results,paste0(opt$sampleid,"_RCTD_results.rds"))
