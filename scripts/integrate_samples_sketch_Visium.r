renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--samplesheet"), type="character", default=NULL, 
              help="A .txt file sample-level meta data. [default= %default]", metavar="character"),
  make_option(c("--geneinfo"), type="character", default="NA", 
              help="Path to gene level annotation file. This is used to add feature level meta data. [default= %default]", metavar="character"),
  make_option(c("--integration_method"), type="character", default="harmony", 
              help="cca, rpca, harmony, fastmnn, or scvi. [default= %default]", metavar="character"),
  make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution parameter used to identify number of clusters. [default= %default]", metavar="double"),
  make_option(c("--cellcycle_correction_flag"), type="character", default="1", 
              help="0 or 1 indicating whether to estimate and correct for cell-cycle effect. If set to 1, also need to specify gene list for S and G2M phase. [default= %default]", metavar="character"),
  make_option(c("--genelist_S_phase"), type="character", default="NA", 
              help="Path to gene list (Gene Symbols) for cell-cycle S-phase, one gene per line. [default= %default]", metavar="character"),
  make_option(c("--genelist_G2M_phase"), type="character", default="NA", 
              help="Path to gene list (Gene Symbols) for cell-cycle G2M-phase, one gene per line. [default= %default]", metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(matrixStats)
library(harmony)
library(BPCells)
library(SeuratObject)
## fix sketchdata verbose issue
library(Seurat)
library(SeuratWrappers)
library(future)
library(doFuture)
library(Banksy)
options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)
registerDoFuture()
plan("multicore", workers = 5)

remove_object <- function(object_name){
  for(i in object_name) assign(i, NULL,envir = .GlobalEnv)
  rm(list = object_name, envir = .GlobalEnv)
  invisible(gc(verbose = FALSE))
}

t1 <- Sys.time()
assay <- "Spatial"
message("perform sample integration using sketch-based method")
seurat_obj_paths <- list.files("./", pattern = "*.rds")
object.list <- foreach(i = 1:length(seurat_obj_paths)) %dopar%{
  object <- readRDS(seurat_obj_paths[i])
  write_matrix_dir(mat = object[[assay]]$counts, dir = paste0("./on_disk_mat_", unique(object$sampleid)))
  mat <- open_matrix_dir(dir = paste0("./on_disk_mat_", unique(object$sampleid)))
  list(metadata = object@meta.data, mat = mat)
}## using normalizedata instead of sctransform
## currently, sctransform has issues with sketch-based analysis
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
integrated_obj <- CreateSeuratObject(counts = counts.mat, min.cells=0,min.features=0, meta.data = metadata)

remove_object(c("counts.mat", "metadata"))

DefaultAssay(integrated_obj) <- assay
integrated_obj <- NormalizeData(integrated_obj)
integrated_obj[[assay]] <- split(integrated_obj[[assay]], f = integrated_obj$sampleid)
# findvariablefeatures and sketchdata are both very slow in speed
# especially sketchdata, stuck on "Calcuating Leverage Score"
# https://github.com/satijalab/seurat/issues/8127
integrated_obj <- FindVariableFeatures(integrated_obj)
if(opt$cellcycle_correction_flag == "1"){
  message('remove cell-cycle related genes from integration features')
  s.features <- read.delim(opt$genelist_S_phase, header = FALSE)$V1
  g2m.features <- read.delim(opt$genelist_G2M_phase, header = FALSE)$V1
  VariableFeatures(integrated_obj) <- setdiff(VariableFeatures(integrated_obj),c(s.features, g2m.features))
}

integrated_obj <- SketchData(object = integrated_obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(integrated_obj) <- "sketch"
integrated_obj <- FindVariableFeatures(integrated_obj, verbose = F)
if(opt$cellcycle_correction_flag == "1"){
  VariableFeatures(integrated_obj) <- setdiff(VariableFeatures(integrated_obj),c(s.features, g2m.features))
}
integrated_obj <- ScaleData(integrated_obj, verbose = F)

integrated_obj <- RunBanksy(integrated_obj, lambda = 0.2, assay = 'sketch', slot = 'data', features = 'variable',group = 'sampleid', split.scale = FALSE, k_geom = 50, dimx = 'x', dimy = 'y')

integrated_obj <- RunPCA(integrated_obj, assay = 'BANKSY',verbose = F, npcs = 30)
DefaultAssay(integrated_obj) <- "sketch"

method_list <- data.frame(name = c("cca", "rpca", "harmony", "fastmnn", "scvi"),
                          function_name = c("CCAIntegration", "RPCAIntegration", "HarmonyIntegration", "FastMNNIntegration", "scVIIntegration"))
integrated_obj <- IntegrateLayers(
  object = integrated_obj, 
  method = get(method_list$function_name[method_list$name == opt$integration_method]),
  orig.reduction = "pca", 
  new.reduction = opt$integration_method,
  verbose = FALSE
)

integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30, reduction = opt$integration_method, assay = "BANKSY")
integrated_obj <- FindClusters(integrated_obj, resolution = opt$resolution, graph.name = "BANKSY_snn", cluster.name = "seurat_clusters_sketch")
integrated_obj <- RunUMAP(integrated_obj, dims = 1:30, reduction = opt$integration_method, assay = "BANKSY")
integrated_obj <- RunTSNE(integrated_obj, dims = 1:30, reduction = opt$integration_method, assay = "BANKSY")

integrated_obj <- ProjectIntegration(object = integrated_obj, 
                                 sketched.assay = "sketch",
                                 assay = assay, 
                                 reduction = opt$integration_method)
integrated_obj <- ProjectData(object = integrated_obj, 
                          sketched.assay = "sketch",
                          assay = assay, 
                          sketched.reduction = opt$integration_method,
                          full.reduction = paste0(opt$integration_method,".full"), dims = 1:30, refdata = list(seurat_clusters = "seurat_clusters_sketch"))
integrated_obj <- RunUMAP(integrated_obj, reduction = paste0(opt$integration_method, ".full"), dims = 1:30, reduction.name = "umap.full",reduction.key = "UMAPfull_")
integrated_obj <- RunTSNE(integrated_obj, reduction = paste0(opt$integration_method, ".full"), dims = 1:30, reduction.name = "tsne.full",reduction.key = "TSNEfull_")
t2 <- Sys.time()
message(paste0('Total time spent for integration: ', as.numeric(t2-t1),' ',units(t2-t1)))

integrated_obj$seurat_clusters <- factor(integrated_obj$seurat_clusters, levels=0:(length(unique(integrated_obj$seurat_clusters))-1))
Idents(integrated_obj) <- integrated_obj$seurat_clusters
if(length(unique(integrated_obj$condition)) > 1) integrated_obj$seurat_clusters_diff <- paste0(integrated_obj$seurat_clusters, "_", integrated_obj$condition)
## gene meta info was removed after integration
## add gene meta info for SCT assay
if(opt$geneinfo != "NA"){
  geneinfo <- data.table::fread(opt$geneinfo, header = TRUE,data.table = FALSE)
  colnames(geneinfo)[1:2] <- c("gene_id", "gene_name")
  geneinfo$gene_name_unique <- make.unique(geneinfo$gene_name)
  rownames(geneinfo) <- geneinfo$gene_name_unique
  integrated_obj@assays[["sketch"]] <- AddMetaData(integrated_obj@assays[["sketch"]], metadata = geneinfo)
}

sampleinfo <- read.delim(opt$samplesheet, header = TRUE)
if(ncol(sampleinfo) > 3){
  additional_metadata <- sampleinfo[match(integrated_obj$sampleid,sampleinfo$sampleid),3:ncol(sampleinfo)]
  rownames(additional_metadata) <- colnames(integrated_obj)
  integrated_obj <- AddMetaData(integrated_obj, metadata = additional_metadata)
}

write.table(data.frame(clusternum=unique(integrated_obj$seurat_clusters)), "seurat_clusters.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(data.frame(condition = integrated_obj$condition, 
                       seurat_clusters = integrated_obj$seurat_clusters), "seurat_clusters_condition.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
integrated_obj[[assay]] <- JoinLayers(integrated_obj[[assay]])
integrated_obj[["sketch"]] <- JoinLayers(integrated_obj[["sketch"]])
DefaultAssay(integrated_obj) <- assay

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

for(i in 1:length(integrated_obj@assays[[assay]]@layers[["counts"]]@matrix@matrix_list))
  integrated_obj@assays[[assay]]@layers[["counts"]]@matrix@matrix_list[[i]] <- update_dir(integrated_obj@assays[[assay]]@layers[["counts"]]@matrix@matrix_list[[i]], "./on_disk_mat")
for(i in 1:length(integrated_obj@assays[[assay]]@layers[["data"]]@matrix@matrix_list))  integrated_obj@assays[[assay]]@layers[["data"]]@matrix@matrix_list[[i]] <- update_dir(integrated_obj@assays[[assay]]@layers[["data"]]@matrix@matrix_list[[i]], "./on_disk_mat")

message("saving integrated object")
saveRDS(integrated_obj, "integrated_obj.rds")