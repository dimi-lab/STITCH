seurat_obj <- readRDS("integrated_obj.rds")

seurat_obj$seurat_clusters <- droplevels(seurat_obj$seurat_clusters)

CellDimPlot(seurat_obj, group.by = "seurat_clusters", reduction = "umap", theme_use = "theme_blank")

CellDimPlot(seurat_obj, group.by = "seurat_clusters", reduction = "umap", theme_use = "theme_blank", split.by = "sampleid")

CellStatPlot(seurat_obj, stat.by = "condition", group.by = "seurat_clusters", stat_type = "count", position = "dodge", label = TRUE)  + 
  labs(x='cluster number', y='number of cells') + theme(axis.text = element_text(size=15), axis.title = element_text(size=18,face="bold"))

CellStatPlot(seurat_obj, stat.by = "condition", group.by = "seurat_clusters", stat_type = "percent", position = "dodge", label = TRUE)  + 
  labs(x='cluster number', y='number of cells') + theme(axis.text = element_text(size=15), axis.title = element_text(size=18,face="bold"))
