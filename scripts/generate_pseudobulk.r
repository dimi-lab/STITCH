renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--seurat_obj"), type="character", default=NULL, 
              help="Path to seurat object file. [default= %default]", metavar="character"),
  make_option(c("--data_type"), type="character", default=NULL, 
              help="data type. [default= %default]", metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(BPCells)
library(Seurat)
options(stringsAsFactors = FALSE)

seurat_obj <- readRDS(opt$seurat_obj)
if(opt$data_type == "scRNAseq") assay <- "RNA" else if(opt$data_type == "Visium") assay <- "Spatial"
seurat_obj_pb <- AggregateExpression(seurat_obj, assays = assay, return.seurat = TRUE, group.by = c("condition", "sampleid", "seurat_clusters"))

write.table(data.frame(condition = seurat_obj_pb$condition, 
                       seurat_clusters = seurat_obj_pb$seurat_clusters), "seurat_clusters_condition.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

saveRDS(seurat_obj_pb, "seurat_obj_pb.rds")