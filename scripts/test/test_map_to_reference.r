# library(BPCells)
# library(Seurat)
# options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)

## test SCT integrated object as reference
opt <- list(
  ## integrated object based on v3 integration for project1
  reference = "testdata/integrated_obj_scRNAseq_V3.rds",
  query = "project2/integrate/integrated_obj.rds",
  reference_reduction = "pca",
  normalization_method = "SCT",
  refdata = "cell_type",
  prediction_assay = "1",
  reduction_model = "umap"
)

reference <- readRDS(opt$reference)
query <- readRDS(opt$query)
ref_vars <- strsplit(opt$refdata, ",")[[1]]
refdata <- as.list(ref_vars)
names(refdata) <- ref_vars
reference$cell_type <- paste0('cluster',reference$seurat_clusters)

anchor <- FindTransferAnchors(
  reference = reference,
  query = query,
  reduction = "pcaproject",
  reference.reduction = "pca",
  normalization.method = "SCT")

query <- TransferData(
  anchorset = anchor, 
  reference = reference,
  query = query,
  refdata = refdata,
  prediction.assay = TRUE
)

query <- IntegrateEmbeddings(
  anchorset = anchor,
  reference = reference,
  query = query, 
  new.reduction.name = paste0("ref.", "pca")
)

query <- ProjectUMAP(
  query = query, 
  query.reduction = paste0("ref.", "pca"), 
  reference = reference, 
  reference.reduction = "pca", 
  reduction.model = "umap"
)

######## test sketch assay and lognormalize
## to make on disk mat accessible

opt <- list(
  reference = "integrate/integrated_obj.rds",
  query = "loader/GEX_CNT_3-19/seurat_obj.rds",
  reference_reduction = "pca",
  normalization_method = "RNA",
  refdata = "cell_type",
  prediction_assay = "1",
  reduction_model = "umap"
)

reference <- readRDS(opt$reference)
query <- readRDS(opt$query)
ref_vars <- strsplit(opt$refdata, ",")[[1]]
refdata <- as.list(ref_vars)
names(refdata) <- ref_vars
reference$cell_type <- paste0('cluster',reference$seurat_clusters)

DefaultAssay(query) <- "RNA"
query <- NormalizeData(query, verbose = FALSE)
query <- FindVariableFeatures(query)
query <- ScaleData(query)
rownames(query) <- make.unique(query@assays$RNA@meta.data$gene_name)

DefaultAssay(reference) <- "sketch"

anchor <- FindTransferAnchors(
  reference = reference,
  query = query,
  reduction = "pcaproject",
  reference.reduction = "pca",
  normalization.method = "LogNormalize")

query <- TransferData(
  anchorset = anchor, 
  reference = reference,
  query = query,
  refdata = refdata,
  prediction.assay = TRUE
)

query <- IntegrateEmbeddings(
  anchorset = anchor,
  reference = reference,
  query = query, 
  new.reduction.name = paste0("ref.", "pca")
)

query <- ProjectUMAP(
  query = query, 
  query.reduction = paste0("ref.", "pca"), 
  reference = reference, 
  reference.reduction = "pca", 
  reduction.model = "umap"
)
