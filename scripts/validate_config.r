library('optparse')
option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL, 
              help="Full path to config file. [default= %default]", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

source("helper_functions.r")
message('parse config file')
config_values <- parse_config(opt$config)
opt <- c(opt, config_values)

var_names <- c("config", "sampleinfo", "species", "output", "data_type", "resolution", "vismethod", "genelist", "controlvar", "experimentvar", "memoryefficient", "allinone", "integration_analysis","covariate_list", "test.use", "QC_only","adaptive_cutoff", "MT_cutoff","HB_cutoff","nFeature_cutoff","nCount_cutoff", "help")
message_list <- c("Start validating provided config file")
if(!setequal(names(opt), var_names))
  message_list <- c(message_list,paste0("variable names do not match predefined variables: ", setdiff(names(opt), var_names), setdiff(var_names, names(opt))))

if(!opt$species %in% c("human", "mouse", "NA"))
  message_list <- c(message_list, 'species has to be one of: human, mouse or NA')

if(! dir.exists(opt$output)){
  message_list <- c(message_list, 'output directory does not exist!')
} else if(! R.utils::isAbsolutePath(opt$output)){
  message_list <- c(message_list, 'output path is not absolute path!')
}

if(!opt$data_type %in% c("scRNAseq", "Visium"))
  message_list <- c(message_list, 'data_type has to be one of: scRNAseq, Visium')

if(!opt$vismethod %in% c('tsne','umap'))
  message_list <- c(message_list, 'Visulization method has to be either tsne or umap!')

if(opt$genelist != "NA" & (!file.exists(opt$genelist)))
  message_list <- c(message_list, 'genelist has to be either NA or a file!')

if(!opt$memoryefficient %in% c('0','1','2','3'))
  message_list <- c(message_list, 'Memory options can only be chosen from 0,1,2,3')

if(!opt$allinone %in% c('0','1'))
  message_list <- c(message_list, 'allinone option can only be chosen from 0,1')

if(!opt$integration_analysis %in% c("0", "1"))
  message_list <- c(message_list, 'integration_analysis option can only be chosen from 0,1')

if(!opt$test.use %in% c("wilcox", "wilcox_limma", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
  message_list <- c(message_list, 'test.use has to be one of wilcox, wilcox_limma, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2')

if((opt$covariate_list != "NA") & (!opt$test.use %in% c("LR", "negbinom", "poisson", "MAST")))
  message_list <- c(message_list, 'covariate_list specified, but test.use is not one of LR, negbinom, poisson or MAST')

if(!opt$QC_only %in% c("0", "1"))
  message_list <- c(message_list, "QC_only has to be one of 0 or 1")

if(!opt$adaptive_cutoff %in% c("0", "1"))
  message_list <- c(message_list, "adaptive_cutoff has to be one of 0 or 1")

for(i in c("resolution", "MT_cutoff", "HB_cutoff", "nFeature_cutoff", "nCount_cutoff")) opt[[i]] <- as.numeric(opt[[i]])

if(!file.exists(opt$sampleinfo)){
  message_list <- c(message_list, 'sampleinfo does not exist!')
  message(paste0(message_list, collapse = "\n"))
  quit()
}

sampleinfo <- read.delim(opt$sampleinfo,header=T,stringsAsFactors=F)
if((opt$controlvar=='NA' | opt$experimentvar=='NA') & length(unique(sampleinfo[,2]))>1)
  message_list <- c(message_list, 'Control or experimental conditions set to 0, but multiple conditions specified in sampleinfo file!')

if(opt$controlvar!='NA' & opt$experimentvar!='NA' & !(opt$controlvar %in% sampleinfo[,2] & opt$experimentvar %in% sampleinfo[,2]))
  message_list <- c(message_list, 'case or control specified, but are not defined in sampleinfo file')

if(opt$covariate_list != "NA"){
  opt$covariate_list <- strsplit(opt$covariate_list, ',')[[1]]
  if(!all(opt$covariate_list %in% colnames(sampleinfo)))
    message_list <- c(message_list, 'covariate_list specified, but are not present as column names of sampleinfo file')
}

if(length(message_list) == 1) message("The provided config file looks good!") else message(paste0(message_list, collapse = "\n"))