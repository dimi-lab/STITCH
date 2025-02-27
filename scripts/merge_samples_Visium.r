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
              help="k_geom parameter for Banksy. Local neighborhood size. Larger values will yield larger domains. [default= %default]", metavar="numeric"),
  make_option(c("--parallel_strategy"), type="character", default=NULL, 
              help="Parallel strategy for future. See help page for plan for details. [default= %default]", metavar="character"),
  make_option(c("--nworkers"), type="integer", default=NULL, 
              help="Number of workers/cpus used for future. [default= %default]", metavar="integer")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(matrixStats)
library(SeuratObject)
## fix sketchdata verbose issue
library(Seurat)
library(SeuratWrappers)
library(future)
library(doFuture)
library(Banksy)
options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)
registerDoFuture()
plan(opt$parallel_strategy, workers = as.integer(opt$nworkers))

remove_object <- function(object_name){
  for(i in object_name) assign(i, NULL,envir = .GlobalEnv)
  rm(list = object_name, envir = .GlobalEnv)
  invisible(gc(verbose = FALSE))
}

t1 <- Sys.time()
message("perform sample merging")
seurat_obj_paths <- list.files("./", pattern = "*.rds")
object.list <- foreach(i = seurat_obj_paths) %dopar% readRDS(i)
#object.features <- unique(unlist(lapply(object.list, function(i) VariableFeatures(i))))
object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
if(opt$cellcycle_correction_flag == "1"){
  message('remove cell-cycle related genes from variable features')
  s.features <- read.delim(opt$genelist_S_phase, header = FALSE)$V1
  g2m.features <- read.delim(opt$genelist_G2M_phase, header = FALSE)$V1
  object.features <- setdiff(object.features,c(s.features, g2m.features))
}

assay <- ifelse(opt$norm_dimreduc == "SCT", "SCT", "Spatial")
if(length(object.list) == 1) merged_obj <- object.list[[1]] else merged_obj <- merge(object.list[[1]], y = object.list[2:length(object.list)], merge.data = TRUE)
DefaultAssay(merged_obj) <- assay
VariableFeatures(merged_obj) <- object.features
remove_object("object.list")
if(opt$norm_dimreduc == "LogNormalize"){
  message("run ScaleData to scale across all samples")
  merged_obj <- ScaleData(merged_obj)
}

if(opt$spatial_cluster == "Banksy"){
  if(assay == "Spatial") merged_obj[[assay]] <- JoinLayers(merged_obj[[assay]])
  merged_obj <- RunBanksy(merged_obj, lambda = opt$lambda, assay = assay, slot = 'data', features = 'variable',group = 'sampleid', split.scale = FALSE, k_geom = opt$k_geom, dimx = 'x', dimy = 'y')
  
  merged_obj <- RunPCA(merged_obj, verbose = FALSE, assay = 'BANKSY', npcs = 30, features = VariableFeatures(merged_obj, assay = assay))
  merged_obj <- FindNeighbors(merged_obj, dims = 1:30, assay = "BANKSY")
  merged_obj <- RunUMAP(merged_obj, dims = 1:30, assay = "BANKSY")
  merged_obj <- RunTSNE(merged_obj, dims = 1:30, assay = "BANKSY")
  merged_obj <- FindClusters(merged_obj, verbose = FALSE, resolution = opt$resolution, graph.name = "BANKSY_snn")
} else {
  merged_obj <- RunPCA(merged_obj, verbose = FALSE, npcs = 30)
  merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
  merged_obj <- RunUMAP(merged_obj, dims = 1:30)
  merged_obj <- RunTSNE(merged_obj, dims = 1:30)
  merged_obj <- FindClusters(merged_obj, verbose = FALSE, resolution = opt$resolution)
}

DefaultAssay(merged_obj) <- assay
t2 <- Sys.time()
message(paste0('Total time spent for merging: ', as.numeric(t2-t1),' ',units(t2-t1)))
Idents(merged_obj) <- factor(Idents(merged_obj),levels=0:(length(unique(merged_obj$seurat_clusters))-1))
if(length(unique(merged_obj$condition)) > 1) merged_obj$seurat_clusters_diff <- paste0(merged_obj$seurat_clusters, "_", merged_obj$condition)
#message('normalizing median UMI counts across SCT models')
#merged_obj <- PrepSCTFindMarkers(merged_obj)
if(opt$geneinfo != "NA"){
  geneinfo <- data.table::fread(opt$geneinfo, header = TRUE,data.table = FALSE)
  colnames(geneinfo)[1:2] <- c("gene_id", "gene_name")
  geneinfo$gene_name_unique <- make.unique(geneinfo$gene_name)
  rownames(geneinfo) <- geneinfo$gene_name_unique
  merged_obj[[assay]] <- AddMetaData(merged_obj[[assay]], metadata = geneinfo)
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
message("saving merged object")
saveRDS(merged_obj, "merged_obj.rds")