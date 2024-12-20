library('optparse')
option_list <- list(                                    
  make_option(c("--seurat_obj"), type="character", default=NULL, 
              help="Path to reference seurat object in .rds format. [default= %default]", metavar="character"),
  make_option(c("--assay"), type="character", default=NULL, 
              help="Assay to use. [default= %default]", metavar="character"),
  make_option(c("--res_low"), type="double", default=0.01, 
              help="Lower end resolution value. [default= %default]", metavar="double"),
  make_option(c("--res_high"), type="double", default=2, 
              help="Higher end resolution value. [default= %default]", metavar="double"),
  make_option(c("--res_step"), type="double", default=0.005, 
              help="Step size for resolution. [default= %default]", metavar="double")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(BPCells)
library(Seurat)
library(ggplot2)
#library(ROGUE)
library(future)
library(doFuture)
registerDoFuture()
## multisession do not work somehow
plan("multicore", workers = 20)
options(stringsAsFactors = F, future.globals.maxSize=500*1024^3)
source("./rogue.r")

## SCT counts processed by PrepSCTFindMarkers does not seem 
## to generate meaningful results
## SCT assay overall do not perform better than RNA count directly
## even using aggregated object
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

rogue_list <- foreach(i = 1:length(res_list)) %dopar% {
  tmp <- FindClusters(object = seurat_obj, resolution = res_list[i])
  ### some times encouter issue with https://github.com/PaulingLiu/ROGUE/issues/10
  ### need to increase span
  ### use rogue_fix to fix issue when no significant gene can be found
  rogue.res <- rogue_fix(expr, labels = tmp$seurat_clusters, samples = tmp$sampleid, platform = "UMI", span = 0.8)
  rogue.res
}

df <- data.frame(n_clusters = n_clusters,
                 resolution = res_list,
                 rogue_score = sapply(rogue_list, function(i) median(unlist(i),na.rm = TRUE)))

gp <- ggplot() + geom_point(aes(x= n_clusters, y = rogue_score), data = df) + theme_classic() + labs(x = "# of clusters", y = "ROGUE score") + theme(axis.title = element_text(face = "bold", size =15), axis.text = element_text(size = 8))
ggsave("./rogue_score.pdf", gp, width = 8, height = 8, useDingbats=FALSE)

# ## automatically detect saturation point
# df$first_diff <- c(NA, diff(df$rogue_score))
# df$second_diff <- c(NA, diff(df$first_diff))
# elbow_point <- df$n_clusters[which.max(abs(df$second_diff))]
# elbow_point
# 
# #devtools::install_github("etam4260/kneedle")
# library(kneedle, lib.loc = "/research/bsi/projects/staff_analysis/m182980/Rlib4.2.2")
# smerc
# kneedle(df$n_clusters, df$rogue_score)
# 
# devtools::install_github("likelet/VSOLassoBag")
# library(VSOLassoBag, lib.loc = "/research/bsi/projects/staff_analysis/m182980/Rlib4.2.2")
# df1 <- data.frame(variable = n_clusters, Frequency = df$rogue_score)
# VSOLassoBag::kneedle(df1)
# 
# library(smerc, lib.loc = "/research/bsi/projects/staff_analysis/m182980/Rlib4.2.2")
# elbow_point(df$n_clusters, df$rogue_score)

saveRDS(list(res_list, n_clusters,rogue_list), "./rogue_objs_aggregated_obj_by_sample.rds")
