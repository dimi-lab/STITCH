renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
  make_option(c("--workflowpath"), type="character", default=NULL, 
              help="workflowpath. [default= %default]", metavar="character"),
  make_option(c("--geneinfo"), type="character", default=NULL, 
              help="geneinfo [default= %default]", metavar="character"),
  make_option(c("--fc"), type="character", default=NULL, 
              help="fc [default= %default]", metavar="character"),
  make_option(c("--pval"), type="character", default=NULL, 
              help="pval [default= %default]", metavar="character"),
  make_option(c("--pval_flag"), type="character", default=NULL, 
              help="pval_flag [default= %default]", metavar="character"),
  make_option(c("--pct"), type="character", default=NULL, 
              help="pct [default= %default]", metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

options(stringsAsFactors = FALSE)
for(i in c("fc","pval","pct")) opt[[i]] <- as.numeric(opt[[i]])

message(paste0("combine tables for marker genes"))
command <- paste('Rscript', file.path(opt$workflowpath,'scripts/combine_tables.r'), 
                 '--geneinfo', opt$geneinfo, 
                 '--reorder_sheets','1',
                 '--fc', opt$fc,
                 '--direction','up',
                 '--pval', opt$pval,
                 '--pval_flag', opt$pval_flag,
                 '--pct', opt$pct/100, 
                 '--prefix', "marker_",sep = " "
)
system(command, wait = TRUE)

if(!file.exists("./NO_FILE")){
  message(paste0("combine tables for diff genes for integration analysis"))
  command <- paste('Rscript', file.path(opt$workflowpath,'scripts/combine_tables.r'),
                   '--geneinfo', opt$geneinfo,
                   '--reorder_sheets','1',
                   '--fc', opt$fc,
                   '--direction','both',
                   '--pval', opt$pval,
                   '--pval_flag', opt$pval_flag,
                   '--pct', opt$pct/100, 
                   '--prefix', "diff_",sep = " "
  )
  system(command, wait = TRUE)
}
