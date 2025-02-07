renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(    
  make_option(c("--samplesheet"), type="character", default=NULL, 
              help="A .txt file sample-level meta data. [default= %default]", metavar="character"),
  make_option(c("--geneinfo"), type="character", default="NA", 
              help="Path to gene level annotation file. This is used to add feature level meta data. [default= %default]", metavar="character"),
  make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution parameter used to identify number of clusters. [default= %default]", metavar="double"),
  make_option(c("--cellcycle_correction_flag"), type="character", default="1", 
              help="0 or 1 indicating whether to estimate and correct for cell-cycle effect. If set to 1, also need to specify gene list for S and G2M phase. [default= %default]", metavar="character"),
  make_option(c("--genelist_S_phase"), type="character", default="NA", 
              help="Path to gene list (Gene Symbols) for cell-cycle S-phase, one gene per line. [default= %default]", metavar="character"),
  make_option(c("--genelist_G2M_phase"), type="character", default="NA", 
              help="Path to gene list (Gene Symbols) for cell-cycle G2M-phase, one gene per line. [default= %default]", metavar="character"),
  make_option(c("--norm_dimreduc"), type="character", default=NULL, 
              help="normalization method for dimension reduction. [default= %default]", metavar="character"),
  make_option(c("--spatial_cluster"), type="character", default=NULL, 
              help="Spatial clustering algorithm, Seurat or Banksy. [default= %default]", metavar="character"),
  make_option(c("--lambda"), type="numeric", default=0.2, 
              help="lambda parameter for Banksy. Influence of the neighborhood. Larger values yield more spatially coherent domains. [default= %default]", metavar="numeric"),
  make_option(c("--k_geom"), type="numeric", default=50, 
              help="k_geom parameter for Banksy. Local neighborhood size. Larger values will yield larger domains. [default= %default]", metavar="numeric")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(matrixStats)
library(BPCells)
library(SeuratObject)
## fix sketchdata verbose issue
library(Seurat)
library(future)
library(doFuture)
options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)
registerDoFuture()
plan("multisession", workers = 5)

remove_object <- function(object_name){
  for(i in object_name) assign(i, NULL,envir = .GlobalEnv)
  rm(list = object_name, envir = .GlobalEnv)
  invisible(gc(verbose = FALSE))
}

t1 <- Sys.time()
assay <- "Spatial"
message("perform sample merging using sketch-based method")
seurat_obj_paths <- list.files("./", pattern = "*.rds")
object.list <- foreach(i = 1:length(seurat_obj_paths)) %dopar%{
  object <- readRDS(seurat_obj_paths[i])
  write_matrix_dir(mat = object[[assay]]$counts, dir = paste0("./on_disk_mat_", unique(object$sampleid)))
  mat <- open_matrix_dir(dir = paste0("./on_disk_mat_", unique(object$sampleid)))
  list(metadata = object@meta.data, mat = mat)
}
## using normalizedata instead of sctransform
## currently, sctransform has issues with sketch-based analysis, https://github.com/satijalab/seurat/issues/7336
## collapse does not work for merge
data.list <- lapply(object.list, function(i) i[[2]])
names(data.list) <- gsub(".rds","",basename(seurat_obj_paths))
metadata.list <- lapply(object.list, function(i) i[[1]])
metadata <- Reduce(rbind, metadata.list)
sample_merge <- CreateSeuratObject(counts = data.list, meta.data = metadata)
sample_merge[[assay]] <- JoinLayers(sample_merge[[assay]])
if(dir.exists("./on_disk_mat")) unlink("./on_disk_mat", recursive = TRUE)
message("create on-disk matrix using merged object")
write_matrix_dir(mat = sample_merge[[assay]]$counts, dir = "./on_disk_mat")
metadata <- sample_merge@meta.data

remove_object(c("sample_merge", "metadata.list", "data.list", "object.list"))

counts.mat <- open_matrix_dir(dir = "./on_disk_mat")
message("create seurat object based on on-disk matrix")
merged_obj <- CreateSeuratObject(counts = counts.mat, min.cells=0,min.features=0, meta.data = metadata)

remove_object(c("counts.mat", "metadata"))

DefaultAssay(merged_obj) <- assay
merged_obj <- NormalizeData(merged_obj)
merged_obj[[assay]] <- split(merged_obj[[assay]], f = merged_obj$sampleid)
merged_obj <- FindVariableFeatures(merged_obj)

if(opt$cellcycle_correction_flag == "1"){
  message('remove cell-cycle related genes from variable features')
  s.features <- read.delim(opt$genelist_S_phase, header = FALSE)$V1
  g2m.features <- read.delim(opt$genelist_G2M_phase, header = FALSE)$V1
  VariableFeatures(merged_obj) <- setdiff(VariableFeatures(merged_obj),c(s.features, g2m.features))
}

merged_obj <- SketchData(object = merged_obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(merged_obj) <- "sketch"
merged_obj <- FindVariableFeatures(merged_obj, verbose = F)
if(opt$cellcycle_correction_flag == "1"){
  VariableFeatures(merged_obj) <- setdiff(VariableFeatures(merged_obj),c(s.features, g2m.features))
}
merged_obj <- ScaleData(merged_obj, verbose = F)

if(opt$spatial_cluster == "Banksy"){
  merged_obj[[assay]] <- JoinLayers(merged_obj[[assay]])
  merged_obj <- RunBanksy(merged_obj, lambda = opt$lambda, assay = 'sketch', slot = 'data', features = 'variable',group = 'sampleid', split.scale = FALSE, k_geom = opt$k_geom, dimx = 'x', dimy = 'y')
  merged_obj <- RunPCA(merged_obj, assay = 'BANKSY',verbose = F, npcs = 30)
  merged_obj <- FindNeighbors(merged_obj, dims = 1:30, reduction = "pca", assay = "BANKSY")
  merged_obj <- FindClusters(merged_obj, resolution = opt$resolution, cluster.name = "seurat_clusters_sketch", graph.name = "BANKSY_snn")
  merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "pca", assay = "BANKSY")
  merged_obj <- RunTSNE(merged_obj, dims = 1:30, reduction = "pca", assay = "BANKSY")
} else {
  merged_obj <- RunPCA(merged_obj, verbose = FALSE, npcs = 30)
  merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
  merged_obj <- RunUMAP(merged_obj, dims = 1:30)
  merged_obj <- RunTSNE(merged_obj, dims = 1:30)
  merged_obj <- FindClusters(merged_obj, verbose = FALSE, resolution = opt$resolution)
}

DefaultAssay(merged_obj) <- assay
merged_obj <- ProjectIntegration(object = merged_obj, 
                                 sketched.assay = "sketch",
                                 assay = assay, 
                                 reduction = "pca")
merged_obj <- ProjectData(object = merged_obj, 
                          sketched.assay = "sketch",
                          assay = assay, 
                          sketched.reduction = "pca",
                          full.reduction = paste0("pca",".full"), dims = 1:30, refdata = list(seurat_clusters = "seurat_clusters_sketch"))
merged_obj <- RunUMAP(merged_obj, reduction = paste0("pca", ".full"), dims = 1:30, reduction.name = "umap.full",reduction.key = "UMAPfull_")
merged_obj <- RunTSNE(merged_obj, reduction = paste0("pca", ".full"), dims = 1:30, reduction.name = "tsne.full",reduction.key = "TSNEfull_")
t2 <- Sys.time()
message(paste0('Total time spent for merging: ', as.numeric(t2-t1),' ',units(t2-t1)))

merged_obj$seurat_clusters <- factor(merged_obj$seurat_clusters, levels=0:(length(unique(merged_obj$seurat_clusters))-1))
Idents(merged_obj) <- merged_obj$seurat_clusters
if(length(unique(merged_obj$condition)) > 1) merged_obj$seurat_clusters_diff <- paste0(merged_obj$seurat_clusters, "_", merged_obj$condition)
## gene meta info was removed after integration
# cannot add gene meta info for on-disk assay
if(opt$geneinfo != "NA"){
  geneinfo <- data.table::fread(opt$geneinfo, header = TRUE,data.table = FALSE)
  colnames(geneinfo)[1:2] <- c("gene_id", "gene_name")
  geneinfo$gene_name_unique <- make.unique(geneinfo$gene_name)
  rownames(geneinfo) <- geneinfo$gene_name_unique
  merged_obj@assays[["sketch"]] <- AddMetaData(merged_obj@assays[["sketch"]], metadata = geneinfo)
}
sampleinfo <- read.delim(opt$samplesheet, header = TRUE)
if(ncol(sampleinfo) > 3){
  additional_metadata <- sampleinfo[match(merged_obj$sampleid,sampleinfo$sampleid),3:ncol(sampleinfo)]
  rownames(additional_metadata) <- colnames(merged_obj)
  merged_obj <- AddMetaData(merged_obj, metadata = additional_metadata)
}
write.table(data.frame(clusternum=unique(merged_obj$seurat_clusters)), "seurat_clusters.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(data.frame(condition = merged_obj$condition, 
                       seurat_clusters = merged_obj$seurat_clusters), "seurat_clusters_condition.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
merged_obj[[assay]] <- JoinLayers(merged_obj[[assay]])
merged_obj[["sketch"]] <- JoinLayers(merged_obj[["sketch"]])
DefaultAssay(merged_obj) <- assay

message('update path to on_disk_mat')
# https://github.com/satijalab/seurat/issues/8691
# https://github.com/satijalab/seurat-object/issues/198
update_dir <- function(object, new_dir) {
  if ("dir" %in% slotNames(object)) {
    slot(object, "dir") <- new_dir
    return(object)
  }
  
  if ("matrix" %in% slotNames(object)) {
    matrix_slot <- slot(object, "matrix")
    updated_matrix_slot <- update_dir(matrix_slot, new_dir)
    slot(object, "matrix") <- updated_matrix_slot
    return(object)
  }
  return(object)
}

for(i in 1:length(merged_obj@assays[[assay]]@layers[["counts"]]@matrix@matrix_list))
  merged_obj@assays[[assay]]@layers[["counts"]]@matrix@matrix_list[[i]] <- update_dir(merged_obj@assays[[assay]]@layers[["counts"]]@matrix@matrix_list[[i]], "./on_disk_mat")
for(i in 1:length(merged_obj@assays[[assay]]@layers[["data"]]@matrix@matrix_list))  merged_obj@assays[[assay]]@layers[["data"]]@matrix@matrix_list[[i]] <- update_dir(merged_obj@assays[[assay]]@layers[["data"]]@matrix@matrix_list[[i]], "./on_disk_mat")

message("saving merged object")
saveRDS(merged_obj, "merged_obj.rds")