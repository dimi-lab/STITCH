renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--workflowpath"), type="character", default=NULL, 
              help="workflow path. [default= %default]", metavar="character"),
  make_option(c("--data_type"), type="character", default=NULL, 
              help="data type. [default= %default]", metavar="character"),
  make_option(c("--ambient_RNA_removal_flag"), type="character", default=NULL, 
              help="ambient RNA removal flag. [default= %default]", metavar="character"),
  make_option(c("--doublet_removal_flag"), type="character", default=NULL, 
              help="doublet removal flag. [default= %default]", metavar="character"),
  make_option(c("--adaptive_cutoff_flag"), type="character", default=NULL, 
              help="adaptive cutoff flag. [default= %default]", metavar="character"),
  make_option(c("--mt_cutoff"), type="character", default=NULL, 
              help="mt cutoff. [default= %default]", metavar="character"),
  make_option(c("--hb_cutoff"), type="character", default=NULL, 
              help="hb cutoff. [default= %default]", metavar="character"),
  make_option(c("--nFeature_cutoff"), type="character", default=NULL, 
              help="nFeature cutoff. [default= %default]", metavar="character"),
  make_option(c("--nCount_cutoff"), type="character", default=NULL, 
              help="nCount cutoff. [default= %default]", metavar="character"),
  make_option(c("--nCell_cutoff"), type="character", default=NULL, 
              help="nCell cutoff. [default= %default]", metavar="character"),
  make_option(c("--geneinfo"), type="character", default=NULL, 
              help="gene info. [default= %default]", metavar="character"),
  make_option(c("--cellcycle_correction_flag"), type="character", default=NULL, 
              help="cellcycle correction flag. [default= %default]", metavar="character"),
  make_option(c("--genelist_S_phase"), type="character", default=NULL, 
              help="genelist S phase. [default= %default]", metavar="character"),
  make_option(c("--genelist_G2M_phase"), type="character", default=NULL, 
              help="genelist G2M phase. [default= %default]", metavar="character"),
  make_option(c("--min_median_umi"), type="character", default=NULL, 
              help="Path to min_median_umi file [default= %default]", metavar="integer"),
  make_option(c("--norm_dimreduc"), type="character", default=NULL, 
              help="normalization method for dimension reduction. [default= %default]", metavar="character"),
  make_option(c("--norm_diff"), type="character", default=NULL, 
              help="normalization method for differential testing. [default= %default]", metavar="character"),
  make_option(c("--sampleid"), type="character", default=NULL, 
              help="sample id. [default= %default]", metavar="integer"),
  make_option(c("--condition"), type="character", default=NULL, 
              help="condition. [default= %default]", metavar="integer"),
  make_option(c("--secondary_output"), type="character", default=NULL, 
              help="Path to secondary output. [default= %default]", metavar="integer")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(matrixStats)
library(SoupX)
library(scDblFinder)
library(ggplot2)
library(Seurat)
library(patchwork)
options(stringsAsFactors = FALSE)

for(i in c("resolution", "mt_cutoff", "hb_cutoff", "nFeature_cutoff", "nCount_cutoff","nCell_cutoff")) opt[[i]] <- as.numeric(opt[[i]])
source(file.path(opt$workflowpath, "scripts/helper_functions.r"))

if(opt$data_type == 'scRNAseq') {
  seurat_obj <- qcsample_scRNAseq(
    sampleid = opt$sampleid,
    condition = opt$condition,
    secondary_output = opt$secondary_output,
    stage = "loader",
    workflowpath = opt$workflowpath,
    ambient_RNA_removal_flag = opt$ambient_RNA_removal_flag,
    doublet_removal_flag = opt$doublet_removal_flag,
    adaptive_cutoff_flag = opt$adaptive_cutoff_flag,
    nCount_cutoff = opt$nCount_cutoff,
    nFeature_cutoff = opt$nFeature_cutoff,
    mt_cutoff = opt$mt_cutoff,
    hb_cutoff = opt$hb_cutoff,
    nCell_cutoff = opt$nCell_cutoff,
    geneinfo = opt$geneinfo,
    cellcycle_correction_flag = opt$cellcycle_correction_flag,
    genelist_S_phase = opt$genelist_S_phase,
    genelist_G2M_phase = opt$genelist_G2M_phase,
    min_median_umi = opt$min_median_umi,
    norm_dimreduc = opt$norm_dimreduc,
    norm_diff = opt$norm_diff
  )
}


if(opt$data_type == 'Visium') {
  seurat_obj <- qcsample_Visium(
    sampleid = opt$sampleid,
    condition = opt$condition,
    secondary_output = opt$secondary_output,
    stage = "loader",
    workflowpath = opt$workflowpath,
    adaptive_cutoff_flag = opt$adaptive_cutoff_flag,
    nCount_cutoff = opt$nCount_cutoff,
    nFeature_cutoff = opt$nFeature_cutoff,
    mt_cutoff = opt$mt_cutoff,
    hb_cutoff = opt$hb_cutoff,
    nCell_cutoff = opt$nCell_cutoff,
    geneinfo = opt$geneinfo,
    cellcycle_correction_flag = opt$cellcycle_correction_flag,
    genelist_S_phase = opt$genelist_S_phase,
    genelist_G2M_phase = opt$genelist_G2M_phase,
    min_median_umi = opt$min_median_umi,
    norm_dimreduc = opt$norm_dimreduc,
    norm_diff = opt$norm_diff
  )
}

message("saving seurat object")
saveRDS(seurat_obj, paste0(opt$sampleid,".rds"))
