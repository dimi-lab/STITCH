renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--workflowpath"), type="character", default=NULL, 
              help="workflow path. [default= %default]", metavar="character"),
  make_option(c("--samplesheet"), type="character", default=NULL, 
              help="sample sheet. [default= %default]", metavar="character"),
  make_option(c("--data_type"), type="character", default=NULL, 
              help="data type. [default= %default]", metavar="character"),
  make_option(c("--resolution"), type="character", default=NULL, 
              help="resolution. [default= %default]", metavar="character"),
  make_option(c("--geneinfo"), type="character", default=NULL, 
              help="gene info. [default= %default]", metavar="character"),
  make_option(c("--cellcycle_correction_flag"), type="character", default=NULL, 
              help="cellcycle correction flag. [default= %default]", metavar="character"),
  make_option(c("--genelist_S_phase"), type="character", default=NULL, 
              help="genelist S phase. [default= %default]", metavar="character"),
  make_option(c("--genelist_G2M_phase"), type="character", default=NULL, 
              help="genelist G2M phase. [default= %default]", metavar="character"),
  make_option(c("--integration_method"), type="character", default=NULL, 
              help="integration method. [default= %default]", metavar="character"),
  make_option(c("--sketch_flag"), type="character", default=NULL, 
              help="sketch flag. [default= %default]", metavar="character"),
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

options(stringsAsFactors = FALSE)

for(i in c("resolution", "lambda", "k_geom")) opt[[i]] <- as.numeric(opt[[i]])

if(opt$data_type == 'scRNAseq') {
  if(opt$sketch_flag == "0"){
    command <- paste("Rscript", file.path(opt$workflowpath,"/scripts/integrate_samples_scRNAseq.r"),
                     "--samplesheet", opt$samplesheet,
                     "--geneinfo", opt$geneinfo,
                     "--integration_method",opt$integration_method,
                     "--resolution", opt$resolution,
                     "--cellcycle_correction_flag", opt$cellcycle_correction_flag,
                     "--genelist_S_phase", opt$genelist_S_phase,
                     "--genelist_G2M_phase", opt$genelist_G2M_phase,
                     "--norm_dimreduc", opt$norm_dimreduc, sep = " "
    )
    system(command, wait = TRUE)
  } else{
    command <- paste("Rscript", file.path(opt$workflowpath,"/scripts/integrate_samples_sketch_scRNAseq.r"),
                     "--samplesheet", opt$samplesheet,
                     "--geneinfo", opt$geneinfo,
                     "--integration_method",opt$integration_method,
                     "--resolution", opt$resolution,
                     "--cellcycle_correction", opt$cellcycle_correction,
                     "--genelist_S_phase", opt$genelist_S_phase,
                     "--genelist_G2M_phase", opt$genelist_G2M_phase, sep = " "
    )
    system(command, wait = TRUE)
  }
}

if(opt$data_type == 'Visium') {
  if(opt$sketch_flag == "0"){
    command <- paste("Rscript", file.path(opt$workflowpath,"/scripts/integrate_samples_Visium.r"),
                     "--samplesheet", opt$samplesheet,
                     "--geneinfo", opt$geneinfo,
                     "--integration_method",opt$integration_method,
                     "--resolution", opt$resolution,
                     "--cellcycle_correction_flag", opt$cellcycle_correction_flag,
                     "--genelist_S_phase", opt$genelist_S_phase,
                     "--genelist_G2M_phase", opt$genelist_G2M_phase,
                     "--norm_dimreduc", opt$norm_dimreduc, 
                     "--spatial_cluster", opt$spatial_cluster, 
                     "--lambda", opt$lambda,
                     "--k_geom", opt$k_geom, sep = " "
    )
    system(command, wait = TRUE)
  } else{
    command <- paste("Rscript", file.path(opt$workflowpath,"/scripts/integrate_samples_sketch_Visium.r"),
                     "--samplesheet", opt$samplesheet,
                     "--geneinfo", opt$geneinfo,
                     "--integration_method",opt$integration_method,
                     "--resolution", opt$resolution,
                     "--cellcycle_correction", opt$cellcycle_correction,
                     "--genelist_S_phase", opt$genelist_S_phase,
                     "--genelist_G2M_phase", opt$genelist_G2M_phase, 
                     "--norm_dimreduc", opt$norm_dimreduc, 
                     "--spatial_cluster", opt$spatial_cluster,
                     "--lambda", opt$lambda,
                     "--k_geom", opt$k_geom, sep = " "
    )
    system(command, wait = TRUE)
  }
}