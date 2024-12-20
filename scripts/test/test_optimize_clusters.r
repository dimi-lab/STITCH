opt <- list(seurat_obj = "project2/integrate/integrated_obj.rds",
            assay = "RNA",
            res_low = 0.01,
            res_high = 2,
            res_step = 0.005)


seurat_obj <- readRDS(opt$seurat_obj)
expr <- as.matrix(seurat_obj[[opt$assay]]$counts)
rownames(expr) <- rownames(seurat_obj[[opt$assay]])
colnames(expr) <- colnames(seurat_obj)
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)

if(length(grep("_snn_res", colnames(seurat_obj@meta.data)))) seurat_obj@meta.data[grep("_snn_res", colnames(seurat_obj@meta.data), value = TRUE)] <- NULL

message("Identify resolutions for unique set of clusters")
res_list <- seq(opt$res_low, opt$res_high, opt$res_step)
tmp <- FindClusters(object = seurat_obj, resolution = res_list, verbose = FALSE)
col_idx <- grep("_snn_res", colnames(tmp@meta.data))[1]
res_list <- res_list[which(!duplicated(apply(tmp@meta.data[,col_idx:ncol(tmp@meta.data)],2,function(i) length(unique(i)))))]
n_clusters <- unique(apply(tmp@meta.data[,col_idx:ncol(tmp@meta.data)],2,function(i) length(unique(i))))
res_list <- res_list[order(n_clusters, decreasing = FALSE)]
n_clusters <- n_clusters[order(n_clusters, decreasing = FALSE)]
res_list <- res_list[n_clusters != 1]
n_clusters <- n_clusters[n_clusters != 1]

