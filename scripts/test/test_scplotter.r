seurat_obj <- readRDS("integrated_obj.rds")

seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj@assays$sketch))

CellDimPlot(seurat_obj,group_by = "seurat_clusters", split_by = "sampleid", reduction = "umap", theme = "theme_blank", highlight =TRUE, label= TRUE, label_insitu = TRUE)

CellDimPlot(seurat_obj_sel,group_by = "seurat_clusters", split_by = "condition", reduction = "umap", theme = "theme_blank", highlight =TRUE)

CellDimPlot(seurat_obj_sel,group_by = "Phase", split_by = "condition", reduction = "umap", theme = "theme_blank", highlight =TRUE)

FeatureStatPlot(seurat_obj, plot_type = "dim", features = c("APOE", "KLF5"), reduction = "umap", theme = "theme_blank", split_by = "condition")

rows_data <- data.frame(genes = sample(rownames(seurat_obj),24),
                        clusters = rep(paste0("cluster",0:7), each =3))


FeatureStatPlot(seurat_obj_heatmap, features = rows_data$genes, ident = "condition", cell_type = "dot", plot_type = "heatmap", name = "Expression Level", dot_size = function(x) sum(x > 0) / length(x), dot_size_name = "Percent Expressed", add_bg = FALSE, rows_data = rows_data, show_row_names = TRUE, rows_split_by = "clusters", cluster_rows = FALSE, cluster_columns = FALSE, row_name_annotation = TRUE, rows_split_palette = "Set2", columns_split_palette = "Set2", columns_split_by = "seurat_clusters", layer = "scale.data")

FeatureStatPlot(seurat_obj_heatmap, features = rows_data$genes, ident = "condition_test", plot_type = "heatmap", name = "Expression Level", rows_data = rows_data, show_row_names = TRUE, rows_split_by = "clusters", cluster_rows = FALSE, row_name_annotation = TRUE, cluster_columns = FALSE, rows_split_palette = "Set2", columns_split_palette = "Set2", layer = "scale.data", columns_split_by = "seurat_clusters")




