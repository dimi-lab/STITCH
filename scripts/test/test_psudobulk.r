library(DESeq2)
library(Seurat)
library(edgeR)
library(limma)
seurat_obj <- readRDS("output_scRNAseq_PBMC_160K/integrate/integrated_obj.rds")

## create artificial condition
seurat_obj$condition_art <- ifelse(seurat_obj$sampleid %in% c("P1","P2","P3","P4"),"condition1", "condition2")
# pseudobulk the counts based on donor-condition-celltype
seurat_obj_pb <- AggregateExpression(seurat_obj, assays = "RNA", return.seurat = TRUE, group.by = c("condition_art", "sampleid", "seurat_clusters"))
seurat_obj_pb$seurat_clusters_condition <- paste(seurat_obj_pb$seurat_clusters, seurat_obj_pb$condition_art, sep = "_")

## compare between conditions
Idents(seurat_obj_pb) <- "seurat_clusters_condition"
bulk_de <- FindMarkers(object = seurat_obj_pb,
                       ident.1 = "3_condition1",
                       ident.2 = "3_condition2",
                       test.use = "DESeq2")
head(bulk_de, n = 15)

## compare between clusters
Idents(seurat_obj_pb) <- "seurat_clusters"
bulk_de_cluster <- FindMarkers(object = seurat_obj_pb,
                       ident.1 = "3",
                       test.use = "DESeq2")
