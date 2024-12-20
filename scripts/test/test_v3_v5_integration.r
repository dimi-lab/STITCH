# ### use cite-seq data
# library(matrixStats)
# library(harmony)
# library(SeuratObject)
# library(Seurat)
# library(SeuratWrappers)
# library(future)
# library(doFuture)
# options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)
# registerDoFuture()
# plan("multisession", workers = 5)

opt <- list(samplesheet = "PBMC_ADT_scRNAseq/sampleinfo.tsv",
            geneinfo = "docs/human_gene_info_2020A.tsv",
            genelist_S_phase = "SPARTA/docs/S_genes_human.tsv",
            genelist_G2M_phase = "docs/G2M_genes_human.tsv",
            resolution = 0.8,
            cellcycle_correction_flag = "1",
            integration_method = "harmony")

reference <- SeuratDisk::LoadH5Seurat("data/multi.h5seurat")

reference[["SCT"]] <- split(reference[["SCT"]], f = reference$donor)

t1 <- Sys.time()
reference <- IntegrateLayers(
  object = reference, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = TRUE,
  features = VariableFeatures(reference)
)
t2 <- Sys.time()
t2 -t1

reference <- JoinLayers(reference)
t1 <- Sys.time()
reference <- RunHarmony(object = reference,
                         assay.use = "SCT",
                         reduction.use = "pca",
                         dims.use = 1:30,
                         group.by.vars = "donor",
                         reduction.save = "rna.harmony.v3",
                         plot_convergence = FALSE)
t2 <- Sys.time()
t2 -t1

integrated_obj <- readRDS("output/integrated_obj.rds")

t1 <- Sys.time()
integrated_obj <- RunHarmony(object = integrated_obj,
                        assay.use = "SCT",
                        reduction.use = "pca",
                        dims.use = 1:30,
                        group.by.vars = "sampleid",
                        reduction.save = "rna.harmony.v3",
                        plot_convergence = FALSE)
t2 <- Sys.time()
t2 -t1
## takes 3 mins

integrated_obj[["SCT"]] <- split(integrated_obj[["SCT"]], f = integrated_obj$sampleid)
DefaultAssay(integrated_obj) <- "SCT"

t1 <- Sys.time()
integrated_obj <- IntegrateLayers(
  object = integrated_obj, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "rna.harmony.v5",
  normalization.method = "SCT",
  verbose = TRUE,
  features = VariableFeatures(integrated_obj)
)
t2 <- Sys.time()
t2 -t1
## get stuck for harmony
## takes 35 mins
## seems to be caused by this issue:
## https://github.com/satijalab/seurat/issues/7879
## reduced to 2 mins after fix

