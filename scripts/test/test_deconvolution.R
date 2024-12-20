# library(spacexr)
# library(Seurat)

## test SCT integrated object as reference
opt <- list(
  ## integrated object based on v3 integration for project1
  reference = "ClusterAnnotated_rename.rds",
  query = "data_merge.rds"
)

reference <- readRDS(opt$reference)
query <- readRDS(opt$query)

reference <- Reference(reference[["SCT"]]$counts, as.factor(reference$ClusterID))

query_sel <- subset(query, sampleid == "S1")
bulk_spatial <- SpatialRNA(GetTissueCoordinates(query_sel,image = "S1"), query_sel@assays$SCT@counts)

### can only use 1 core to avoid parallel bug
myRCTD <- create.RCTD(bulk_spatial, reference, max_cores = 1, CELL_MIN_INSTANCE=20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

saveRDS(myRCTD@results$weights, paste0(outputdir, "/weights_", samplename,".rds"))
