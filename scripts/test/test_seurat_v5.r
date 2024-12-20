###############################
######## test integration pipeline with sctransform
######## without BPcells or sketch pipeline
sample1 <- readRDS("GEX_CNT_3-19/seurat_obj.rds")
sample2 <- readRDS("loader/GEX_IL-17_3-19/seurat_obj.rds")
# sample_list <- lapply(1:20, function(i) {
#   sample1$sampleid <- paste0("sample",i)
#   sample1})
sample_list <- list(sample1, sample2)
integrated_obj <- merge(sample_list[[1]], sample_list[2:length(sample_list)])
VariableFeatures(integrated_obj) <- SelectIntegrationFeatures(sample_list)
integrated_obj <- RunPCA(integrated_obj, npcs = 30, verbose = F)

# integrated_obj <- SCTransform(integrated_obj, method="glmGamPoi", assay = "RNA", variable.features.n = 3000,vst.flavor = "v2",verbose = FALSE,return.only.var.genes = FALSE)
# integrated_obj <- RunPCA(integrated_obj, npcs = 30, verbose = F)
integrated_obj <- IntegrateLayers(
  object = integrated_obj,
  method = FastMNNIntegration,
  normalization.method = "SCT",
  verbose = F,
  orig.reduction = "pca"
)
integrated_obj <- RunUMAP(integrated_obj, dims = 1:30, reduction = "integrated.dr")
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30, reduction = "integrated.dr")
integrated_obj <- FindClusters(integrated_obj, resolution = 2)
gp <- DimPlot(integrated_obj, reduction = "umap", group.by = "seurat_clusters", split.by = "sampleid")
ggsave("integration_test_umap.png", gp, width = 15, height=10)


#datapaths <- paste0(list.files("/research/bsi/archive/PI/Beyder_Arthur_m045616/secondary/s306116.spatialGI/mrnaseq/230427-AB/" ,pattern = "^[A-D]",full.names = TRUE), "/outs/filtered_feature_bc_matrix.h5")

## something wrong with igraph installation
## keeps getting errors of libicui18n.so.58: cannot open shared object file: No such file or directory
#library(igraph, lib.loc = '/usr/local/biotools/rpackages/R-4.2.2-2023-02-01')

# library(matrixStats)
# library(Seurat)
# library(BPCells)
# #library(future)
# library(patchwork)
# #registerDoFuture()
# #plan("multisession", workers = 10)
# options(future.globals.maxSize = 500*1024^3)

sample1 <- readRDS("output_scRNAseq/loader/GEX_CNT_3-19/seurat_obj.rds")
#sample2 <- readRDS("/research/bsi/projects/staff_analysis/m182980/script/seurat/5.1/testdata/091324/output_scRNAseq/loader/GEX_IL-17_3-19//seurat_obj.rds")

### convert a seurat object to on-disk matrix for BPcells
### collapse not functioning now
#sample_merge <- Reduce(merge, lapply(1:20, function(i) sample1))
sample_list <- lapply(1:20, function(i) {
  sample1$sampleid <- paste0("sample",i)
  sample1})
sample_merge <- merge(sample_list[[1]], sample_list[2:length(sample_list)])
rm(list = "sample_list")
gc()
sample_merge[["RNA"]] <- JoinLayers(sample_merge[["RNA"]])
meta.data <- sample_merge@meta.data
write_matrix_dir(mat = sample_merge[["RNA"]]$counts, dir = '100324/on_disk_mat/')
rm(list = "sample_merge")
gc()
counts.mat <- open_matrix_dir(dir = '100324/on_disk_mat/')
seurat_obj <- CreateSeuratObject(counts = counts.mat, min.cells=0,min.features=0, meta.data = meta.data)
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$sampleid)
# format(object.size(seurat_obj), units = "Mb")
# seurat_obj1 <- seurat_obj
# seurat_obj1[["RNA"]]$`counts.GEX_CNT_3-19` <- as(object = seurat_obj[["RNA"]]$`counts.GEX_CNT_3-19`, Class = "dgCMatrix")
# seurat_obj1[["RNA"]]$`counts.GEX_IL-17_3-19` <- as(object = seurat_obj[["RNA"]]$`counts.GEX_IL-17_3-19`, Class = "dgCMatrix")
# format(object.size(seurat_obj1), units = "Mb")

t1 <- Sys.time()
seurat_obj <- SCTransform(seurat_obj, method="glmGamPoi", variable.features.n = 3000,vst.flavor = "v2",verbose = FALSE,return.only.var.genes = FALSE)
t2 <- Sys.time()
t2-t1

seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = F)
seurat_obj <- IntegrateLayers(
  object = seurat_obj,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = FALSE
)
seurat_obj <- IntegrateLayers(
  object = seurat_obj,
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  method = CCAIntegration,
  normalization.method = "SCT",
  verbose = FALSE
)
seurat_obj <- IntegrateLayers(
  object = seurat_obj, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = FALSE
)
seurat_obj <- IntegrateLayers(
  object = seurat_obj, 
  method = FastMNNIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
# obj <- IntegrateLayers(
#   object = obj, method = scVIIntegration,
#   new.reduction = "integrated.scvi",
#   conda_env = "../miniconda3/envs/scvi-env", verbose = FALSE
# )

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "integrated.dr")
seurat_obj <- FindClusters(seurat_obj, resolution = 2)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "integrated.dr")
DimPlot(seurat_obj, group.by = "seurat_clusters", split.by = "sampleid", reduction = "umap")

saveRDS(seurat_obj, "seurat_obj.rds")

############## sketch based integration
#### works for NormalizaData, not work for sctransform
## https://github.com/satijalab/seurat/issues/8428
## 
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)

seurat_obj <- SketchData(
  object = seurat_obj,
  ncells = 2000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch to analyzing the sketched dataset (in-memory)
DefaultAssay(seurat_obj) <- "sketch"
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = F)
seurat_obj <- IntegrateLayers(seurat_obj, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca", dims = 1:30, k.anchor = 20, reference = which(Layers(seurat_obj, search = "data") %in% c("data.GEX_CNT_3-19")),verbose = F)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 2)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, return.model = T)
DimPlot(seurat_obj, label = T, label.size = 3, reduction = "umap", split.by = "sampleid") + NoLegend()

seurat_obj <- ProjectIntegration(object = seurat_obj, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")
seurat_obj <- ProjectData(object = seurat_obj, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(cluster_full = "seurat_clusters"))

seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
                  reduction.key = "UMAP_full_")

seurat_obj <- ProjectData(
  object = seurat_obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:30,
  refdata = list(cluster_full = "seurat_clusters")
)
# now that we have projected the full dataset, switch back to analyzing all cells
DefaultAssay(seurat_obj) <- "RNA"
DimPlot(seurat_obj, label = T, label.size = 3, reduction = "full.umap", group.by = "cluster_full", alpha = 1, split.by = "sampleid") + NoLegend()

# visualize gene expression on the sketched cells (fast) and the full dataset (slower)
DefaultAssay(seurat_obj) <- "sketch"
x1 <- FeaturePlot(seurat_obj, "ENSG00000238009")
DefaultAssay(seurat_obj) <- "RNA"
x2 <- FeaturePlot(seurat_obj, "ENSG00000238009")
x1 | x2

#####################
#### test sketch-analysis for 1 million PBMC cells
parse.data <- open_matrix_anndata_hdf5(
  "PBMC_1M/PBMC_1M.h5ad"
)
write_matrix_dir(mat = parse.data, dir = "PBMC_1M/parse_1m_pbmc_112924")

parse.mat <- open_matrix_dir(dir = "parse_1m_pbmc_112924")
# need to move
metadata <- read.csv("PBMC_1M/cell_metadata_1M_PBMC.csv", header = TRUE)
rownames(metadata) <- as.character(0:(nrow(metadata)-1))
object <- CreateSeuratObject(counts = parse.mat, meta.data = metadata)
object <- NormalizeData(object)
# split assay into 24 layers
object[["RNA"]] <- split(object[["RNA"]], f = object$sample)
# findvariablefeatures and sketchdata are both very slow in speed
# especially sketchdata, stuck on "Calcuating Leverage Score"
# https://github.com/satijalab/seurat/issues/8127
# findvariablefeatures works on layer by layer
# have the variable features saved for each sample during loader step?
# will this break the logic of BPcells?
# not enabled to use multi-process
# perform sketch per-sample, then merge?
# foreach does not work
# node stack overflow
#object_sel <- subset(object, sample %in% c("D_505", "D_506"))
#object_sel <- FindVariableFeatures(object_sel, verbose = TRUE)

t1 <- Sys.time()
object <- FindVariableFeatures(object, verbose = TRUE)
t2 <- Sys.time()
t2 - t1
## 18 mins

## do this instead
## no need to save hvfinfo
## take union of features
t1 <- Sys.time()
layer <- Layers(object = object, search = "counts")
vf_list <- foreach(i = seq_along(layer), .combine = 'c') %dopar% {
  features <- Features(x = object, layer = layer[i])
  list(features[FindVariableFeatures(LayerData(object = object, layer = layer[i], fast = TRUE), verbose = TRUE)$variable])
}
t2 <- Sys.time()
t2 - t1
## 1.5 mins
VariableFeatures(object) <- Reduce(union, vf_list)

## sketchdata is also performed layer by layer
t1 <- Sys.time()
object <- SketchData(object = object, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch", verbose = TRUE)
t2 <- Sys.time()
t2 - t1
# takes 1 hour

object
DefaultAssay(object) <- "sketch"
## this step might be optinal, as SketchData already assigned variable features for each layer
#object <- FindVariableFeatures(object, verbose = TRUE)

t1 <- Sys.time()
object <- ScaleData(object, verbose = TRUE)
t2 <- Sys.time()
t2 - t1
## 2.6 mins

t1 <- Sys.time()
object <- RunPCA(object, verbose = TRUE)
t2 <- Sys.time()
t2 - t1
## 2.3 mins

# integrate the datasets
# object <- IntegrateLayers(object, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
#                           dims = 1:30, k.anchor = 20, reference = which(Layers(object, search = "data") %in% c("data.H_3060")),
#                           verbose = F)

t1 <- Sys.time()
object <- IntegrateLayers(
  object = object, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = TRUE
)
t2 <- Sys.time()
t2 - t1
## 30 mins
## reduced to 2 mins after fix

# cluster the integrated data
object <- FindNeighbors(object, reduction = "harmony", dims = 1:30)
object <- FindClusters(object, resolution = 2)
object <- RunUMAP(object, reduction = "harmony", dims = 1:30, return.model = T, verbose = TRUE)

# you can now rejoin the layers in the sketched assay this is required to perform differential
# expression
#object[["sketch"]] <- JoinLayers(object[["sketch"]])
#c10_markers <- FindMarkers(object = object, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
#head(c10_markers)

# You can now annotate clusters using marker genes.  We performed this step, and include the
# results in the 'sketch.celltype' metadata column

plot.s1 <- DimPlot(object, group.by = "sample", reduction = "umap")
plot.s2 <- DimPlot(object, group.by = "seurat_clusters", reduction = "umap")
plot.s1 + plot.s2 + plot_layout(ncol = 1)

# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
#object[["sketch"]] <- split(object[["sketch"]], f = object$sample)
## seems to perform layer by layer
t1 <- Sys.time()
object <- ProjectIntegration(object = object, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
t2 <- Sys.time()
t2 - t1
## 19 mins

t1 <- Sys.time()
object <- ProjectData(object = object, 
                      sketched.assay = "sketch",
                      assay = "RNA", 
                      sketched.reduction = "harmony",
                      full.reduction = "harmony.full", dims = 1:30, refdata = list(seurat_clusters_full = "seurat_clusters"), verbose = TRUE)
t2 <- Sys.time()
t2 - t1
## 7mins

object <- RunUMAP(object, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",reduction.key = "UMAPfull_")

saveRDS(object, "/object.rds")

object <- readRDS("/object.rds")

p1 <- DimPlot(object, reduction = "umap.full", group.by = "sample", alpha = 0.1)
p2 <- DimPlot(object, reduction = "umap.full", group.by = "seurat_clusters_full", alpha = 0.1)
p1 + p2 + plot_layout(ncol = 1)

## differential expression based on psudo-bulk analysis
bulk <- AggregateExpression(object, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("celltype.full",
                                                                                                     "sample", "disease"))

# each sample is an individual-specific celltype-specific pseudobulk profile
tail(Cells(bulk))

cd14.bulk <- subset(bulk, celltype.full == "CD14 Mono")
Idents(cd14.bulk) <- "disease"
de_markers <- FindMarkers(cd14.bulk, ident.1 = "D", ident.2 = "H", slot = "counts", test.use = "DESeq2",
                          verbose = F)
de_markers$gene <- rownames(de_markers)
ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)


########################################################
########### test sketch-workflow for merging-based analysis
t1 <- Sys.time()
parse.mat <- open_matrix_dir(dir = "../../testdata/100324/PBMC_1M/parse_1m_pbmc")
# need to move
metadata <- read.csv("../../testdata/100324/PBMC_1M/cell_metadata_1M_PBMC.csv", header = TRUE)
rownames(metadata) <- as.character(0:(nrow(metadata)-1))
object <- CreateSeuratObject(counts = parse.mat, meta.data = metadata)
format(object.size(object), units = "Mb")
object <- NormalizeData(object)
# split assay into 24 layers
object[["RNA"]] <- split(object[["RNA"]], f = object$sample)
# findvariablefeatures and sketchdata are both very slow in speed
# especially sketchdata, stuck on "Calcuating Leverage Score"
# https://github.com/satijalab/seurat/issues/8127
object <- FindVariableFeatures(object, verbose = FALSE)
object <- SketchData(object = object, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
object
DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object, verbose = F)
object <- ScaleData(object, verbose = F)
object <- RunPCA(object, verbose = F)

# cluster the integrated data
object <- FindNeighbors(object, dims = 1:30)
object <- FindClusters(object, resolution = 2)
object <- RunUMAP(object, dims = 1:30, return.model = T, verbose = F)
object <- ProjectIntegration(object = object, sketched.assay = "sketch", assay = "RNA", reduction = "pca")
object <- ProjectData(object = object, 
                      sketched.assay = "sketch",
                      assay = "RNA", 
                      sketched.reduction = "pca",
                      full.reduction = "pca.full", dims = 1:30, refdata = list(seurat_clusters_full = "seurat_clusters"))
object <- RunUMAP(object, reduction = "pca.full", dims = 1:30, reduction.name = "umap.full",reduction.key = "UMAPfull_")
t2 <- Sys.time()
message(paste0('Total time spent for merging-sketch: ', as.numeric(t2-t1),' ',units(t2-t1)))
gc()
## takes 2 hours and 9.8G
gp <- DimPlot(object, reduction = "umap", group.by = "sample")
ggsave("dimplot_1Mcells_merge.png", gp, width=30, height =20)

########### without sketch
t1 <- Sys.time()
parse.mat <- open_matrix_dir(dir = "../../testdata/100324/PBMC_1M/parse_1m_pbmc")
# need to move
metadata <- read.csv("../../testdata/100324/PBMC_1M/cell_metadata_1M_PBMC.csv", header = TRUE)
rownames(metadata) <- as.character(0:(nrow(metadata)-1))
object <- CreateSeuratObject(counts = parse.mat, meta.data = metadata)
format(object.size(object), units = "Mb")
object <- NormalizeData(object)
# split assay into 24 layers
object[["RNA"]] <- split(object[["RNA"]], f = object$sample)
# findvariablefeatures and sketchdata are both very slow in speed
# especially sketchdata, stuck on "Calcuating Leverage Score"
# https://github.com/satijalab/seurat/issues/8127
object <- FindVariableFeatures(object, verbose = FALSE)
object <- ScaleData(object, verbose = F)
object <- RunPCA(object, verbose = F)

# cluster the integrated data
object <- FindNeighbors(object, dims = 1:30)
object <- FindClusters(object, resolution = 2)
object <- RunUMAP(object, dims = 1:30, return.model = T, verbose = F)
t2 <- Sys.time()
message(paste0('Total time spent for merging: ', as.numeric(t2-t1),' ',units(t2-t1)))
gc()
## takes 8.6 hours and 62G


