###### PBMC data from https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub
###### around 160K cells
# library(future)
# library(doFuture)
# options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)
# registerDoFuture()
# plan("multicore", workers = 10)
# library("DropletUtils")

library(Seurat)
library(patchwork)

reference <- SeuratDisk::LoadH5Seurat("PBMC_ADT_scRNAseq/data/multi.h5seurat", assays = "SCT")

data.seurat.list <- Seurat::SplitObject(reference, split.by = "donor")
sample.names <- unique(reference$donor)

demultiplex_convert_to_10x <- function(obj, samples, assay, outputdir) {
  if(class(obj[[1]]) != "Seurat") {
    message("WARNING: this rds file does not contain a Seurat object! STOP RUNNING THIS SCRIPT")
    message("Check the data type by running:")
    message("class(obj[[1]])")
    stop()
  }
  if(!dir.exists(file.path(outputdir, "demultiplexed"))) {
    dir.create(file.path(outputdir, "demultiplexed"))
  } else {
    print("WARNING! A demultiplexed directory already exists")
    return()
  }
  foreach(i = 1:length(samples)) %dopar% {
    print(paste0("Converting sample ", samples[i]))
    obj.sub <- obj[[samples[i]]]
    DropletUtils::write10xCounts(path = paste0(outputdir,"/demultiplexed/",samples[i]), x = obj.sub[[assay]]@counts, type = "sparse", version="3")
  }
}

demultiplex_convert_to_10x(obj = data.seurat.list, samples = sample.names, assay = "SCT", outputdir = "PBMC_160K/")

sampleinfo <- data.frame(sampleid = sort(unique(reference$donor)),
                         condition = "PBMC",
                         secondary_output = file.path("PBMC_160K/demultiplexed", sort(unique(reference$donor))))
write.table(sampleinfo, "testdata/sampleinfo_scRNAseq_PBMC_160K.txt", row.names = F, col.names = T, sep = "\t", quote = F)


sampleinfo <- data.frame(sampleid = list.files("SPARTA_test/DLPFC/"),
                         condition = "DLPFC",
                         secondary_output = file.path("SPARTA_test/DLPFC", list.files("/SPARTA_test/DLPFC/")))
write.table(sampleinfo, "testdata/sampleinfo_Visium_DLPFC.txt", row.names = F, col.names = T, sep = "\t", quote = F)

######### integrate vs merge for DLPFC
integrated_obj <- readRDS("output_Visium_DLPFC_SCT/integrate/integrated_obj.rds")
p1 <- SpatialDimPlot(integrated_obj, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151673")
p2 <- SpatialDimPlot(integrated_obj, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151507")

merged_obj <- readRDS("output_Visium_DLPFC_SCT/merge/merged_obj.rds")
p3 <- SpatialDimPlot(merged_obj, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151673")
p4 <- SpatialDimPlot(merged_obj, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151507")

########### lognorm
integrated_obj_lognorm <- readRDS("output_Visium_DLPFC_lognorm//integrate/integrated_obj.rds")
p5 <- SpatialDimPlot(integrated_obj_lognorm, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151673")
p6 <- SpatialDimPlot(integrated_obj_lognorm, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151507")


merge_obj_lognorm <- readRDS("output_Visium_DLPFC_lognorm/merge/merged_obj.rds")
p7 <- SpatialDimPlot(merge_obj_lognorm, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151673")
p8 <- SpatialDimPlot(merge_obj_lognorm, group.by = "seurat_clusters", label = T, repel = T, label.size = 4, pt.size.factor = 2, images = "X151507")

wrap_plots(list(p1,p2,p3,p4,p5,p6,p7,p8),ncol = 2)


