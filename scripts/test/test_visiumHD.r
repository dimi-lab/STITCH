library(Seurat)
options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)

data_dir <- "/research/bsi/projects/PI/tertiary/Ordog_Tamas_m038357/s312975.multi_omics/visiumHD_test/small_intestine/output/"
seurat_obj <- Load10X_Spatial(data.dir = data_dir, slice = "SI", bin.size = 8)
assay <- Assays(seurat_obj)
sampleinfo <- data.frame(sampleid = "SI", condition = "test", secondary_output = data_dir)
metadata <- sampleinfo[rep(1,ncol(seurat_obj)),]
rownames(metadata) <- colnames(seurat_obj)
seurat_obj <- AddMetaData(seurat_obj, metadata)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-|^mt-")
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^HB.*|^hb.*")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, assay = assay)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj)

s.features <- read.delim("docs/S_genes_mouse.tsv", header = FALSE)$V1
g2m.features <- read.delim("docs/G2M_genes_mouse.tsv", header = FALSE)$V1
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.features, g2m.features = g2m.features, set.ident = FALSE)
vars.to.regress <- c('percent.mt','Phase')
seurat_obj <- subset(seurat_obj, nCount_Spatial.008um > 50)

### very slow
t1 <- Sys.time()
seurat_obj <- SCTransform(seurat_obj, method="glmGamPoi", assay = assay, variable.features.n = 3000,vst.flavor = "v2",verbose = FALSE,return.only.var.genes = FALSE)
t2 <- Sys.time()
t2-t1

## this takes a huge amounts of memory
## try sketch the data first
t1 <- Sys.time()
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = assay, features = VariableFeatures(seurat_obj), selection.method = "moransi")
t2 <- Sys.time()
t2-t1

moranfast(x, c1, c2, alternative = "two.sided")

