renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--workflowpath"), type="character", default=NULL, 
              help="workflow path. [default= %default]", metavar="character"),
  make_option(c("--data_type"), type="character", default=NULL, 
              help="data type. [default= %default]", metavar="character"),
  make_option(c("--pseudobulk_flag"), type="character", default=NULL, 
              help="0 or 1 indicating whether to perform pseudo-bulk based analysis. [default= %default]", metavar="character"),
  make_option(c("--control_var"), type="character", default=NULL, 
              help="control_var. [default= %default]", metavar="character"),
  make_option(c("--case_var"), type="character", default=NULL, 
              help="case_var. [default= %default]", metavar="character"),
  make_option(c("--covariate_list"), type="character", default=NULL, 
              help="covariate_list. [default= %default]", metavar="character"),
  make_option(c("--test"), type="character", default=NULL, 
              help="test. [default= %default]", metavar="character"),
  make_option(c("--sketch_flag"), type="character", default=NULL, 
              help="sketch_flag. [default= %default]", metavar="character"),
  make_option(c("--norm_diff"), type="character", default=NULL, 
              help="normalization method for differential testing. [default= %default]", metavar="character"),
  make_option(c("--seurat_obj"), type="character", default=NULL, 
              help="Path to seurat object file. [default= %default]", metavar="character"),
  make_option(c("--seurat_clusters_condition"), type="character", default=NULL, 
              help="A .txt file with two columns, condition and seurat_clusters. [default= %default]", metavar="character"),
  make_option(c("--clusternum"), type="integer", default=NULL, 
              help="cluster number. [default= %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

options(stringsAsFactors = FALSE)

seurat_clusters_condition <- read.delim(opt$seurat_clusters_condition, header = TRUE)
condition <- seurat_clusters_condition$condition
seurat_clusters <- seurat_clusters_condition$seurat_clusters

if(opt$sketch_flag == "1" & opt$pseudobulk_flag == "0"){
  assay <- "sketch"
} else {
  if(opt$norm_diff == "SCT") assay <- "SCT" else{
    if(opt$data_type == "scRNAseq") assay <- "RNA"
    if(opt$data_type == "Visium") assay <- "Spatial"
  }
  }

selectedclusters <- as.numeric(which(!apply(table(condition, seurat_clusters),2,function(i) any(i < 3)))-1)
if(length(unique(condition))==1 | (!opt$clusternum %in% selectedclusters)){
  message(paste0("perform marker gene identification for cluster ", opt$clusternum))
  command <- paste('Rscript', file.path(opt$workflowpath,'scripts/find_markers.r'),
                   '--rdsfile',opt$seurat_obj,
                   '--assay', assay,
                   '--idents', 'seurat_clusters',
                   '--ident.1', opt$clusternum,
                   '--test.use', opt$test,
                   '--outputfile', paste0('./marker_gene_cluster', opt$clusternum,'.txt'), sep = ' '
  )
  system(command, wait = TRUE)
} else {
  message(paste0("perform conserved marker gene identification for cluster ", opt$clusternum))
  command <- paste('Rscript', file.path(opt$workflowpath,'scripts/find_markers.r'),
                   '--rdsfile',opt$seurat_obj,
                   '--assay', assay,
                   '--idents', 'seurat_clusters',
                   '--ident.1', opt$clusternum,
                   '--group.var', 'condition',
                   '--test.use', opt$test,
                   '--outputfile', paste0('./marker_gene_cluster', opt$clusternum,'.txt'), sep = " ")
  system(command, wait = TRUE)}
if(length(unique(condition)) > 1 & (opt$control_var != "NA") & (opt$case_var != "NA")){
  idx <- which(condition %in% c(opt$control_var, opt$case_var))
  selectedclusters_diff <- as.numeric(which(!apply(table(condition[idx], seurat_clusters[idx]),2,function(i) any(i==0)))-1)
  if(opt$clusternum %in% selectedclusters_diff){
    message(paste0("perform differential gene identification for cluster ", opt$clusternum, ' between ', opt$case_var, ' and ', opt$control_var))
    command <- paste('Rscript', file.path(opt$workflowpath,'scripts/find_markers.r'),
                     '--rdsfile',opt$seurat_obj,
                     '--assay', assay,
                     '--idents', 'seurat_clusters_diff',
                     '--ident.1', paste0(opt$clusternum, "_", opt$case_var),
                     '--ident.2', paste0(opt$clusternum,"_",opt$control_var),
                     '--test.use', opt$test,
                     '--covariate.list',opt$covariate_list,
                     '--outputfile', paste0('./diff_gene_cluster', opt$clusternum,'.txt'), sep = " "
    )
    system(command, wait = TRUE)
  }
}