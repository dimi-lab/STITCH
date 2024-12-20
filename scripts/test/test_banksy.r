# library(SeuratWrappers)
# library(Banksy)
# library(Seurat)
# library(ggplot2)

object <- readRDS("work/bc/2198bbab618453fb87ca4b2f105ea2/C1_005.rds")
# set lambda to 0.2 for visium data
object <- RunBanksy(object,
                    lambda = 0.2, verbose = TRUE,
                    assay = "SCT", slot = "data", features = "variable",
                    k_geom = 50
)

object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = VariableFeatures(object, assay = "SCT"), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)

Idents(object) <- "banksy_cluster"
p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4, pt.size.factor = 2)
ggsave("./banksy_spatial_plot.png",p)

banksy_cells <- CellsByIdentities(object)
p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
ggsave("./banksy_spatial_plot_highlight.png",p)


# DefaultAssay(object) <- "SCT"
# object <- RunPCA(object, assay = "SCT", reduction.name = "pca", features = rownames(object), npcs = 30)
# object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
# object <- FindClusters(object, cluster.name = "seurat_cluster", resolution = 0.5)
# 
# Idents(object) <- "seurat_cluster"
# p <- SpatialDimPlot(object, group.by = "seurat_cluster", label = T, repel = T, label.size = 4, pt.size.factor = 2)
# ggsave("./seurat_spatial_plot.png",p)
# clustering is much cleaner than seurat-clusters

# integration
object.list <- foreach(i = c("work/bc/2198bbab618453fb87ca4b2f105ea2/C1_005.rds", "work/bc/2198bbab618453fb87ca4b2f105ea2/D1_005.rds")) %dopar% readRDS(i)
metadata1 <- read.csv("C1-005/outs/spatial/tissue_positions.csv", header = TRUE, row.names = 1)
metadata2 <- read.csv("D1-005/outs/spatial/tissue_positions.csv", header = TRUE, row.names = 1)
object.list[[1]] <- AddMetaData(object.list[[1]], metadata = metadata1)
object.list[[2]] <- AddMetaData(object.list[[2]], metadata = metadata2)

message('select integration features')
object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)

integrated_obj <- merge(object.list[[1]], y = object.list[2:length(object.list)],  merge.data = TRUE)
VariableFeatures(integrated_obj) <- SelectIntegrationFeatures(object.list)

integrated_obj <- RunBanksy(integrated_obj, lambda = 0.2, assay = 'SCT', slot = 'data', features = 'variable',group = 'sampleid', split.scale = FALSE, k_geom = 50, dimx = 'pxl_col_in_fullres', dimy = 'pxl_row_in_fullres')

DefaultAssay(integrated_obj) <- "SCT"
integrated_obj <- RunPCA(integrated_obj, assay = 'BANKSY', npcs = 30, features = VariableFeatures(integrated_obj, assay = "SCT"))
integrated_obj <- RunHarmony(integrated_obj, group.by.vars='sampleid', assay = "BANKSY")
integrated_obj <- RunUMAP(integrated_obj, dims = 1:10, reduction = 'harmony', assay = "BANKSY")
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:10, reduction = 'harmony', assay = "BANKSY")
integrated_obj <- FindClusters(integrated_obj, resolution = 0.5, graph.name = "BANKSY_snn")
p <- SpatialDimPlot(integrated_obj, group.by = "seurat_clusters", pt.size.factor = 2, crop = FALSE)
ggsave("banksy_spatial_plot_integrated.png",p)


##### V5 integration workflow
integrated_obj <- merge(object.list[[1]], y = object.list[2:length(object.list)],  merge.data = TRUE)
VariableFeatures(integrated_obj) <- SelectIntegrationFeatures(object.list)

integrated_obj <- RunBanksy(integrated_obj, lambda = 0.2, assay = 'SCT', slot = 'data', features = 'variable',group = 'sampleid', split.scale = FALSE, k_geom = 50, dimx = 'pxl_col_in_fullres', dimy = 'pxl_row_in_fullres')

integrated_obj <- RunPCA(integrated_obj, assay = 'BANKSY', npcs = 30, features = VariableFeatures(integrated_obj, assay = "SCT"))

DefaultAssay(integrated_obj) <- "SCT"
integrated_obj <- IntegrateLayers(
  object = integrated_obj, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = FALSE,
)

######## evaluate on DLPFC sample 
object <- readRDS("work/22/17739ababc853b15228300062e8501/151673.rds")
object <- RunBanksy(object,
                    lambda = 0.2, verbose = TRUE,
                    assay = "Spatial", slot = "data", features = "variable",
                    k_geom = 50, dimx = "x", dimy = "y"
)

object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = VariableFeatures(object, assay = "Spatial"), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.7)

Idents(object) <- "banksy_cluster"
SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4, pt.size.factor = 2)
### lognormalize seems to generate more meaningful clusters compared to SCT

## load all samples and perform intergration
setwd("work/3d/199e08fcfbebe94efb275ad369b781/")
seurat_obj_paths <- setdiff(list.files("./", pattern = "*.rds"),"integrated_obj.rds")
object.list <- foreach(i = seurat_obj_paths[c(1,11)]) %dopar% readRDS(i)
message('select integration features')
object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)

integrated_obj <- merge(object.list[[1]], y = object.list[2:length(object.list)],  merge.data = TRUE)
VariableFeatures(integrated_obj) <- object.features
assay <- "Spatial"
DefaultAssay(integrated_obj) <- assay
# need to join layer, especially scale.data layer
# for calculating PCA
integrated_obj[[assay]] <- JoinLayers(integrated_obj[[assay]], layers = c("counts", "data", "scale.data"))
integrated_obj <- RunPCA(integrated_obj, verbose = TRUE)
integrated_obj[[assay]] <- split(integrated_obj[[assay]], f = integrated_obj$sampleid)
DefaultAssay(integrated_obj) <- assay
integrated_obj <- IntegrateLayers(
  object = integrated_obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  normalization.method = "LogNormalize",
  verbose = TRUE,
  features = VariableFeatures(integrated_obj, assay = assay)
)
reduction <- "harmony"
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30, reduction = reduction)
integrated_obj <- FindClusters(integrated_obj, resolution = 0.3)
integrated_obj <- RunUMAP(integrated_obj, dims = 1:30, reduction = reduction, return.model = TRUE)
integrated_obj <- RunTSNE(integrated_obj, dims = 1:30, reduction = reduction)
DimPlot(integrated_obj, group.by = "seurat_clusters", split.by = "sampleid")
SpatialDimPlot(integrated_obj, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2)

##### scale data across all samples
integrated_obj1 <- merge(object.list[[1]], y = object.list[2:length(object.list)],  merge.data = TRUE)
VariableFeatures(integrated_obj1) <- object.features
assay <- "Spatial"
DefaultAssay(integrated_obj1) <- assay
integrated_obj1 <- ScaleData(integrated_obj1)
integrated_obj1 <- RunPCA(integrated_obj1, verbose = TRUE)
DefaultAssay(integrated_obj1) <- assay
integrated_obj1 <- IntegrateLayers(
  object = integrated_obj1,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  normalization.method = "LogNormalize",
  verbose = TRUE,
  features = VariableFeatures(integrated_obj1, assay = assay)
)
reduction <- "harmony"
integrated_obj1 <- FindNeighbors(integrated_obj1, dims = 1:30, reduction = reduction)
integrated_obj1 <- FindClusters(integrated_obj1, resolution = 0.3)
integrated_obj1 <- RunUMAP(integrated_obj1, dims = 1:30, reduction = reduction, return.model = TRUE)
integrated_obj1 <- RunTSNE(integrated_obj1, dims = 1:30, reduction = reduction)
DimPlot(integrated_obj1, group.by = "seurat_clusters", split.by = "sampleid")
SpatialDimPlot(integrated_obj1, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2)

#### not much difference for integration analysis
#### for merge-based analysis, scale across all samples actually performs better
#### decide to rerun scaledata after merging








integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30, reduction = "harmony", assay = "BANKSY")
integrated_obj <- FindClusters(integrated_obj, resolution = 0.3, graph.name = "BANKSY_snn")
integrated_obj <- RunUMAP(integrated_obj, dims = 1:30, reduction = "harmony", assay = "BANKSY")

SpatialDimPlot(integrated_obj, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151673")




