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
library(harmony)
library(SeuratObject)
## fix sketchdata verbose issue
library(Seurat)
library(SeuratWrappers)
library(future)
library(doFuture)
### need to use R version 4.4.1
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
seurat_obj_paths <- list.files("./", pattern = "*.rds")
object.list <- foreach(i = seurat_obj_paths) %dopar% readRDS(i)
message('select integration features')
object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
if(opt$cellcycle_correction_flag == "1"){
  message('remove cell-cycle related genes from integration features')
  s.features <- read.delim(opt$genelist_S_phase, header = FALSE)$V1
  g2m.features <- read.delim(opt$genelist_G2M_phase, header = FALSE)$V1
  object.features <- setdiff(object.features,c(s.features, g2m.features))
}

assay <- ifelse(opt$norm_dimreduc == "SCT", "SCT", grep("Spatial", Assays(object.list[[1]]), value = TRUE))
integrated_obj <- merge(object.list[[1]], y = object.list[2:length(object.list)],  merge.data = TRUE)
DefaultAssay(integrated_obj) <- assay
remove_object("object.list")
VariableFeatures(integrated_obj) <- object.features
if(opt$norm_dimreduc == "LogNormalize"){
  message("run ScaleData to scale across all samples")
  integrated_obj <- ScaleData(integrated_obj)
}

method_list <- data.frame(name = c("cca", "rpca", "harmony", "fastmnn", "scvi"),
                          function_name = c("CCAIntegration", "RPCAIntegration", "HarmonyIntegration", "FastMNNIntegration", "scVIIntegration"))
if(opt$spatial_cluster == "Banksy"){
  if(length(grep("Spatial", assay))) integrated_obj[[assay]] <- JoinLayers(integrated_obj[[assay]])
  integrated_obj <- RunBanksy(integrated_obj, lambda = opt$lambda, assay = assay, slot = 'data', features = 'variable',group = 'sampleid', split.scale = FALSE, k_geom = opt$k_geom, dimx = 'x', dimy = 'y')
  integrated_obj <- RunPCA(integrated_obj, assay = 'BANKSY', npcs = 30, features = VariableFeatures(integrated_obj, assay = assay))
  integrated_obj[[assay]] <- split(integrated_obj[[assay]], f = integrated_obj$sampleid)
  DefaultAssay(integrated_obj) <- assay
  integrated_obj <- IntegrateLayers(
    object = integrated_obj, 
    method = get(method_list$function_name[method_list$name == opt$integration_method]),
    orig.reduction = "pca", 
    new.reduction = opt$integration_method,
    normalization.method = opt$norm_dimreduc,
    verbose = TRUE,
    features = VariableFeatures(integrated_obj, assay = assay)
  )
  integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30, reduction = opt$integration_method, assay = "BANKSY")
  integrated_obj <- FindClusters(integrated_obj, resolution = opt$resolution, graph.name = "BANKSY_snn")
  integrated_obj <- RunUMAP(integrated_obj, dims = 1:30, reduction = opt$integration_method, assay = "BANKSY")
  integrated_obj <- RunTSNE(integrated_obj, dims = 1:30, reduction = opt$integration_method, assay = "BANKSY")
} else {
  integrated_obj <- RunPCA(integrated_obj, npcs = 30, verbose = F)
  integrated_obj <- IntegrateLayers(
  object = integrated_obj, 
  method = get(method_list$function_name[method_list$name == opt$integration_method]),
  orig.reduction = "pca", 
  new.reduction = opt$integration_method,
  normalization.method = opt$norm_dimreduc,
  verbose = TRUE,
  features = VariableFeatures(integrated_obj, assay = assay)
)
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30, reduction = opt$integration_method)
integrated_obj <- FindClusters(integrated_obj, resolution = opt$resolution)
integrated_obj <- RunUMAP(integrated_obj, dims = 1:30, reduction = opt$integration_method, return.model = TRUE)
integrated_obj <- RunTSNE(integrated_obj, dims = 1:30, reduction = opt$integration_method)
}

t2 <- Sys.time()
message(paste0('Total time spent for integration: ', as.numeric(t2-t1),' ',units(t2-t1)))

Idents(integrated_obj) <- factor(Idents(integrated_obj),levels=0:(length(unique(integrated_obj$seurat_clusters))-1))
if(length(unique(integrated_obj$condition)) > 1) integrated_obj$seurat_clusters_diff <- paste0(integrated_obj$seurat_clusters, "_", integrated_obj$condition)
#message('normalizing median UMI counts across SCT models')
#integrated_obj <- PrepSCTFindMarkers(integrated_obj)
## gene meta info was removed after integration
## add gene meta info for SCT assay
if(opt$geneinfo != "NA"){
  geneinfo <- data.table::fread(opt$geneinfo, header = TRUE,data.table = FALSE)
  colnames(geneinfo)[1:2] <- c("gene_id", "gene_name")
  geneinfo$gene_name_unique <- make.unique(geneinfo$gene_name)
  rownames(geneinfo) <- geneinfo$gene_name_unique
  integrated_obj@assays[[assay]] <- AddMetaData(integrated_obj@assays[[assay]], metadata = geneinfo)
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
message("saving integrated object")
saveRDS(integrated_obj, "integrated_obj.rds")