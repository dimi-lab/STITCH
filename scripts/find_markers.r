renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--rdsfile"), type="character", default=NULL, 
              help="Full path to seurat object in .rds format. [default= %default]", metavar="character"),
  make_option(c("--assay"), type="character", default="SCT", 
              help="assay to use. [default= %default]", metavar="character"),
  make_option(c("--idents"), type="character", default=NULL, 
              help="Character value to set Idents of the seurat object. [default= %default]", metavar="character"),
  make_option(c("--ident.1"), type="character", default=NULL, 
              help="Character value specifying group1. [default= %default]", metavar="character"), 
  make_option(c("--ident.2"), type="character", default=NULL, 
              help="Character value specifying group2. If NULL, use all other cells for comparison. [default= %default]", metavar="character"),
  make_option(c("--group.var"), type="character", default=NULL, 
              help="Character value specifying the group variable. [default= %default]", metavar="character"),
  make_option(c("--test.use"), type="character", default=NULL, 
              help="Character value specifying method for statistical test. See help page for FindMarkers for details. [default= %default]", metavar="character"),
  make_option(c("--covariate.list"), type="character", default=NULL, 
              help="Characters seperated by ',' indicating covariates that need to be adjusted for differential analysis, e.g., age,gender. [default= %default]", metavar="character"),
  make_option(c("--outputfile"), type="character", default=NULL, 
              help="Path to output file in .txt format. [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(matrixStats)
library(BPCells)
library(Seurat)
#library(future)
#library(doFuture)
#options(future.globals.maxSize=100*1024^3,stringsAsFactors = FALSE)
#registerDoFuture()
options(stringsAsFactors = FALSE)
#message(paste0("Number of total workers: ", parallelly::availableCores()))
#plan("multicore", workers = 10)

message('loading rds file')
seurat_obj <- readRDS(opt$rdsfile)
seurat_obj <- SetIdent(seurat_obj, value = opt$idents)
if(!"ident.2" %in% names(opt)) ident.2 <- NULL else ident.2 <- opt$ident.2
if((!"covariate.list" %in% names(opt)) || (opt$covariate.list == "NA")) covariate.list <- NULL else covariate.list <- strsplit(opt$covariate.list,",")[[1]]

seurat_obj[[opt$assay]] <- JoinLayers(seurat_obj[[opt$assay]])

if(!"group.var" %in% names(opt)){
  clusterout <- FindMarkers(seurat_obj, ident.1 = opt$ident.1, ident.2 = ident.2, verbose = TRUE,min.cells.group=1,min.pct=0,min.cells.feature=0,logfc.threshold=0,assay = opt$assay,slot='data', latent.vars = covariate.list, test.use = opt$test.use)
} else {
  ## use default function for now
  clusterout <- FindConservedMarkers(seurat_obj, ident.1 = opt$ident.1, grouping.var = opt$group.var, verbose = TRUE,min.cells.group=1,min.pct=0,min.cells.feature=0,logfc.threshold=0,assay = opt$assay,slot = 'data',latent.vars = covariate.list, test.use = opt$test.use)
  #clusterout <- FindConservedMarkers_fix(seurat_obj, ident.1 = opt$ident.1, grouping.var = opt$group.var, verbose = TRUE,min.cells.group=1,min.pct=0,min.cells.feature=0,logfc.threshold=0,assay = 'SCT',slot = 'data',latent.vars = covariate.list, test.use = opt$test.use)
}
write.table(clusterout,file=opt$outputfile,sep='\t',row.names=T,col.names=T,quote=FALSE)