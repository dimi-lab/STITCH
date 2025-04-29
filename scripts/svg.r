renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--sampleid"), type="character", default=NULL, 
              help="sample id. [default= %default]", metavar="integer"),
  make_option(c("--norm_dimreduc"), type="character", default=NULL, 
              help="normalization method for dimension reduction. [default= %default]", metavar="character"),
  make_option(c("--svg_method"), type="character", default=NULL, 
              help="Method for selection of SVGs, either moransi or markvariogram. [default= %default]", metavar="character")
    )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(Rfast2)
library(Seurat)

seurat_obj <- readRDS(paste0(opt$sampleid,".rds"))
if(opt$norm_dimreduc == "SCT") assay <- "SCT" else assay <- grep("Spatial",Assays(seurat_obj), value = TRUE)
DefaultAssay(seurat_obj) <- assay
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = assay, features = VariableFeatures(seurat_obj), selection.method = opt$svg_method)
if(assay == "SCT"){
  col_sel <- grep(opt$svg_method, colnames(seurat_obj[[assay]]@meta.features), ignore.case = TRUE, value = TRUE)
  svg_res <- seurat_obj[[assay]]@meta.features[rownames(seurat_obj[[assay]]@meta.features) %in% VariableFeatures(seurat_obj),col_sel]
  rownames(svg_res) <- rownames(seurat_obj[[assay]]@meta.features)[rownames(seurat_obj[[assay]]@meta.features) %in% VariableFeatures(seurat_obj)]
} else if(assay == "Spatial"){
  col_sel <- grep(opt$svg_method, colnames(seurat_obj[[assay]]@meta.data), ignore.case = TRUE, value = TRUE)
  svg_res <- seurat_obj[[assay]]@meta.data[rownames(seurat_obj[[assay]]@features) %in% VariableFeatures(seurat_obj),col_sel]
  rownames(svg_res) <- rownames(seurat_obj[[assay]]@features)[rownames(seurat_obj[[assay]]@features) %in% VariableFeatures(seurat_obj)]
}

write.table(svg_res, paste0(opt$sampleid,"_svg_results.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)