library(Rfast2)
library(Seurat)

seurat_obj <- Load10X_Spatial("../../../../../SPARTA_test/DLPFC/151507/")
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)

seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "SCT", features = VariableFeatures(seurat_obj), selection.method = "moransi")
svg <- seurat_obj[["SCT"]]@meta.features[rownames(seurat_obj[["SCT"]]@meta.features) %in% VariableFeatures(seurat_obj),]

DefaultAssay(seurat_obj) <- "Spatial"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, features = VariableFeatures(seurat_obj), selection.method = "moransi")
