library('optparse')
option_list <- list(
  make_option(c("--query"), type="character", default=NULL, 
              help="Path to query seurat object. [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(BPCells)
library(Seurat)
#library(spacexr)
options(future.globals.maxSize = 500*1024^3,stringsAsFactors = FALSE)

query <- readRDS(opt$query)
weightmat <- do.call(rbind, lapply(list.files("./", "weights_.*.rds"), function(i) readRDS(i)))
weightmat <- weightmat[match(colnames(query),rownames(weightmat)),]
query$deconvolution_cell_type <- colnames(weightmat)[apply(weightmat, 1, which.max)]
query[["deconvolution_predictions"]] <- CreateAssay5Object(data = t(weightmat))
saveRDS(query,gsub(".rds","_deconvolution.rds",basename(opt$query)))