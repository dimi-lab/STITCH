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
              help="Path to gene list (Ensembl IDs) for cell-cycle S-phase, one gene per line. [default= %default]", metavar="character"),
  make_option(c("--genelist_G2M_phase"), type="character", default="NA", 
              help="Path to gene list (Ensembl IDs) for cell-cycle G2M-phase, one gene per line. [default= %default]", metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(matrixStats)
library(Seurat)
library(future)
library(doFuture)
library(harmony)
options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)
registerDoFuture()
plan("multisession", workers = 5)

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
if(opt$integration_method =="harmony"){
  message('Generating merged data for harmony')
  integrated_obj <- merge(object.list[[1]], y = object.list[2:length(object.list)],  project = "seurat_harmony", merge.data = TRUE)
  rm(list = c("object.list"))
  gc()
  VariableFeatures(integrated_obj) <- object.features
  integrated_obj <- RunPCA(object = integrated_obj, assay = "SCT",verbose = FALSE ,npcs = 30)
  message('Running harmony')
  integrated_obj <- RunHarmony(object = integrated_obj,
                               assay.use = "SCT",
                               reduction.use = "pca",
                               dims.use = 1:30,
                               group.by.vars = "sampleid",
                               plot_convergence = FALSE)
  integrated_obj <- RunTSNE(integrated_obj,assay='SCT',reduction = "harmony", dims = 1:30)
  integrated_obj <- RunUMAP(object = integrated_obj, assay = "SCT", reduction = "harmony", dims = 1:30)
  integrated_obj <- FindNeighbors(object = integrated_obj, assay = "SCT", reduction = "harmony", dims = 1:30)
  integrated_obj <- FindClusters(object = integrated_obj, resolution = opt$resolution)
} else if(opt$integration_method %in% c("cca", "rpca")){
  message('prepare for SCT integration')
  object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = object.features, verbose = FALSE,assay = 'SCT')
  message('identify integrate anchors')
  if(opt$integration_method=="rpca") {
    object.list <- lapply(X = object.list, FUN = RunPCA, verbose = FALSE, features = object.features)
    mergeAnchor <- FindIntegrationAnchors(object.list = object.list, dims = 1:30,normalization.method = "SCT",anchor.features = object.features,reference = which.max(sapply(object.list,function(i) ncol(i))),reduction = "rpca")
  } else mergeAnchor <- FindIntegrationAnchors(object.list = object.list, dims = 1:30,normalization.method = 'SCT',anchor.features = object.features,verbose = FALSE, reduction = "cca")
  rm(list = c("object.list"))
  gc()
  message('perform integration')
  ### seems to cause issues for v5.1
  ### Warning: Layer counts isn't present in the assay object; returning NULL
  integrated_obj <- IntegrateData(anchorset = mergeAnchor,normalization.method = "SCT",verbose = FALSE,dims = 1:30)
  rm(list = c("mergeAnchor"))
  gc()
  DefaultAssay(integrated_obj) <- "integrated"
  # Run the standard workflow for visualization and clustering
  integrated_obj <- RunPCA(integrated_obj, npcs = 30, verbose = FALSE, features = VariableFeatures(integrated_obj))
  # t-SNE and Clustering
  integrated_obj <- RunTSNE(integrated_obj, reduction = "pca", dims = 1:30)
  integrated_obj <- RunUMAP(integrated_obj,reduction = "pca",dims=1:30)
  integrated_obj <- FindNeighbors(integrated_obj, reduction = "pca", dims = 1:30)
  integrated_obj <- FindClusters(integrated_obj, resolution = opt$resolution)
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
  integrated_obj@assays$SCT <- AddMetaData(integrated_obj@assays$SCT, metadata = geneinfo)
}
sampleinfo <- read.delim(opt$samplesheet, header = TRUE)
if(ncol(sampleinfo) > 3){
  additional_metadata <- sampleinfo[match(integrated_obj$sampleid,sampleinfo$sampleid),3:ncol(sampleinfo)]
  rownames(additional_metadata) <- colnames(integrated_obj)
  integrated_obj <- AddMetaData(integrated_obj, metadata = additional_metadata)
}
write.table(data.frame(clusternum=levels(integrated_obj$seurat_clusters)), "seurat_clusters.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(data.frame(condition = integrated_obj$condition, 
                       seurat_clusters = integrated_obj$seurat_clusters), "seurat_clusters_condition.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
message("saving integrated object")
saveRDS(integrated_obj, "integrated_obj.rds")