# library(spacexr)
# library(Seurat)

## test SCT integrated object as reference
opt <- list(
  ## integrated object based on v3 integration for project1
  reference = "ClusterAnnotated_rename.rds",
  query = "data_merge.rds"
)

reference <- readRDS(opt$reference)
query <- readRDS(opt$query)

reference <- Reference(reference[["SCT"]]$counts, as.factor(reference$ClusterID))

query_sel <- subset(query, sampleid == "S1")
bulk_spatial <- SpatialRNA(GetTissueCoordinates(query_sel,image = "S1"), query_sel@assays$SCT@counts)

### can only use 1 core to avoid parallel bug
myRCTD <- create.RCTD(bulk_spatial, reference, max_cores = 1, CELL_MIN_INSTANCE=20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

saveRDS(myRCTD@results$weights, paste0(outputdir, "/weights_", samplename,".rds"))

## snRNA-seq DLPFC data
## singlecellexperiment object downloaded from https://research.libd.org/globus/
## DLPFC_snRNAseq
library(Seurat)
library(HDF5Array)
se_obj <- loadHDF5SummarizedExperiment(dir="../../../../../SPARTA_test/DLPFC_snRNAseq/", prefix="")
countdata <- as.matrix(assay(se_obj, "counts"))
logcountsdata <- as.matrix(assay(se_obj, "logcounts"))
rownames(countdata) <- make.unique(rownames(countdata))
rownames(logcountsdata) <- make.unique(rownames(logcountsdata))
seurat_obj <- CreateSeuratObject(counts = countdata, meta.data = as.data.frame(colData(se_obj)), min.cells = 0, min.features = 0)
LayerData(seurat_obj[["RNA"]], layer = "data") <- logcountsdata
seurat_obj@reductions[['PCA']] <- CreateDimReducObject(embeddings = reducedDim(se_obj, "GLMPCA_approx"), key = "PC_", assay = "RNA")
seurat_obj@reductions[['HARMONY']] <- CreateDimReducObject(embeddings = reducedDim(se_obj, "HARMONY"), key = "HARMONY_", assay = "RNA")
seurat_obj@reductions[['UMAP']] <- CreateDimReducObject(embeddings = reducedDim(se_obj, "UMAP"), key = "UMAP_", assay = "RNA")

DimPlot(seurat_obj, reduction = "UMAP", group.by = "cellType_broad_hc", split.by = "BrNum", ncol = 5)
DimPlot(seurat_obj, reduction = "UMAP", group.by = "layer_annotation", split.by = "BrNum", ncol = 5)

## recalculate PCA, harmony and umap
## as feature loading is not present in the singlecellexperiment class
## and it's required for transferanchors function
## umap model is not saved either
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$BrNum)
seurat_obj <- seurat_obj %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

seurat_obj <- IntegrateLayers(
  object = seurat_obj, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  normalization.method = "LogNormalize",
  verbose = TRUE
)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", return.model = TRUE, reduction.name = "umap", dims = 1:30)

DimPlot(seurat_obj, reduction = "umap", group.by = "cellType_broad_hc", split.by = "BrNum", ncol = 5)
DimPlot(seurat_obj, reduction = "umap", group.by = "layer_annotation", split.by = "BrNum", ncol = 5)

seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])

saveRDS(seurat_obj, "testdata/process_dir_032825/DLPFC_snRNAseq.rds")
