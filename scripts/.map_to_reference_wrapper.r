library('optparse')
option_list <- list(
  make_option(c("--config"), type="character", default=NULL, 
              help="Path to config file. [default= %default]", metavar="character"),
  make_option(c("--query"), type="character", default=NULL, 
              help="Path to query seurat object. [default= %default]", metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

parse_config <- function(config){
  config <- readLines(config,warn=F)
  config <- grep('^#|^ ',config,value = T,invert = T)
  config <- config[config!=""]
  config_values <- lapply(strsplit(config,'='), function(i) i[2])
  names(config_values) <- sapply(strsplit(config,'='), function(i) i[1])
  return(config_values)
}

options(stringsAsFactors = FALSE)
message('parse config file')
config_values <- parse_config(opt$config)
opt <- c(opt, config_values)
for(i in c("resolution", "MT_cutoff", "HB_cutoff", "nFeature_cutoff", "nCount_cutoff","nCell_cutoff", "pct")) opt[[i]] <- as.numeric(opt[[i]])

if(opt$sketch_flag == "1") assay <- "sketch" 
if(opt$normalization_method_dimreduc == "SCT") assay <- "SCT" else if(opt$normalization_method_dimreduc == "LogNormalize") {
  if(opt$data_type == "scRNAseq") assay <- "RNA" else assay <- "Spatial"
}

command <- paste('Rscript', file.path(opt$workflowpath,'scripts/map_to_reference.r'), 
                 '--reference', opt$reference_scRNAseq, 
                 '--query',opt$query,
                 '--reference_assay',opt$reference_assay,
                 '--query_assay',assay,
                 '--reference_reduction',opt$reference_reduction,
                 '--normalization_method', opt$normalization_method_dimreduc,
                 '--refdata', opt$refdata,
                 '--prediction_assay', "1",
                 '--reduction_model',"umap", sep = " "
)
system(command, wait = TRUE)