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
  make_option(c("--reference_reduction"), type="character", default=NULL, 
              help="Name of dimensional reduction to use from the reference. [default= %default]", metavar="character"),
  make_option(c("--normalization_method"), type="character", default=NULL, 
              help="Name of normalization method used: LogNormalize or SCT. [default= %default]", metavar="character"),
  make_option(c("--refdata"), type="character", default=NULL, 
              help="Character or characters seperated by ',' indicating reference data to use for mapping. Must exist in the metadata of the reference seurat object. [default= %default]", metavar="character"),
  make_option(c("--prediction_assay"), type="character", default=NULL, 
              help="0 or 1 indicating whether to add prediction assay to the seurat object. [default= %default]", metavar="character"),
  make_option(c("--reduction_model"), type="character", default=NULL, 
              help="DimReduc object to use from the reference data. [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(BPCells)
library(Seurat)
options(stringsAsFactors = FALSE)

reference <- readRDS(opt$reference)
query <- readRDS(opt$query)
DefaultAssay(reference) <- opt$reference_assay
DefaultAssay(query) <- opt$query_assay
ref_vars <- strsplit(opt$refdata, ",")[[1]]
refdata <- as.list(ref_vars)
names(refdata) <- ref_vars

message("Performing mapping to reference")

### error when setting reference.reduction to harmony
### https://github.com/satijalab/seurat/issues/8132
### SCT assay for integrated reference is not compatible with 
### FindTransferAnchors https://github.com/satijalab/seurat/issues/9353
### need to run V3 integration workflow with integratedata to create
### integrated assay, use method 'rpca'

anchor <- FindTransferAnchors(
  reference = reference,
  query = query,
  reduction = "pcaproject",
  reference.reduction = opt$reference_reduction,
  normalization.method = opt$normalization_method)

query <- TransferData(
  anchorset = anchors, 
  reference = reference,
  query = query,
  refdata = refdata,
  prediction.assay = as.logical(as.numeric(opt$prediction_assay)))

query <- IntegrateEmbeddings(
  anchorset = anchor,
  reference = reference,
  query = query, 
  new.reduction.name = paste0("ref.", opt$reference_reduction)
)

query <- ProjectUMAP(
  query = query, 
  query.reduction = paste0("ref.", opt$reference_reduction), 
  reference = reference, 
  reference.reduction = opt$reference_reduction, 
  reduction.model = opt$reduction_model
)

message("saving seurat object")
saveRDS(query,gsub(".rds","_mapped_to_ref.rds",basename(opt$query)))