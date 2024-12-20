library('optparse')
option_list <- list(                                    
  make_option(c("--reference"), type="character", default=NULL, 
              help="Path to reference seurat object in .rds format. [default= %default]", metavar="character"),
  make_option(c("--query"), type="character", default=NULL, 
              help="Path to query seurat object. [default= %default]", metavar="character"),
  make_option(c("--reference_assay"), type="character", default=NULL, 
              help="Assay to use to perform mapping for reference. [default= %default]", metavar="character"),
  make_option(c("--query_assay"), type="character", default=NULL, 
              help="Assay to use to perform mapping for query [default= %default]", metavar="character"),
  make_option(c("--refdata"), type="character", default=NULL, 
              help="Character or characters seperated by ',' indicating reference data to use for mapping. Must exist in the metadata of the reference seurat object. [default= %default]", metavar="character"),
  make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample image name for deconvolution. [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(BPCells)
library(Seurat)
#library(spacexr)
options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)

reference <- readRDS(opt$reference)
query <- readRDS(opt$query)

reference <- Reference(reference[[opt$reference_assay]]$counts, as.factor(reference@meta.data[[opt$refdata]]))

query_sel <- subset(query, sampleid == opt$sampleid)
bulk_spatial <- SpatialRNA(GetTissueCoordinates(query_sel,image = opt$sampleid), query_sel[[opt$query_assay]]$counts)

### can only use 1 core to avoid parallel bug
myRCTD <- create.RCTD(bulk_spatial, reference, max_cores = 1, CELL_MIN_INSTANCE=min(table(reference@cell_types)))
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

saveRDS(myRCTD@results$weights, paste0("weights_", opt$sampleid,".rds"))