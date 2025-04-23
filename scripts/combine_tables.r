renv::load(Sys.getenv("PROJECT_DIR"))
library('optparse')
option_list <- list(
    make_option(c("--geneinfo"), type="character", default=NULL, 
              help="Full path to gene level information file to be added to each cluster level table. The first column should be primary gene identifier (usually ensembl ids). The second column should be secondary gene identifier (usually gene symbols). [default= %default]", metavar="character"),
  make_option(c("--reorder_sheets"), type="character", default=NULL, 
              help="1 or 0 to indicate whether the tables should be reordered based on the numeric value after 'cluster'. Set to 1 to enable reorder the sheet names when combining to a .xlsx table. [default= %default]", metavar="character"),
  make_option(c("--fc"), type="double", default=NULL, 
              help="Fold change cutoff to select DEGs. [default= %default]", metavar="double"),
  make_option(c("--direction"), type="character", default="both", 
              help="Direction of fold change to select DEGs. Can be either up (up-regulation) or both (both up and down-regulation). [default= %default]", metavar="character"),
  make_option(c("--pval"), type="double", default=NULL, 
              help="P value cutoff to select DEGs. [default= %default]", metavar="double"),
  make_option(c("--pval_flag"), type="character", default=1, 
              help="0 or 1 to indicate whether to use Bonferroni adjusted p value. [default= %default]", metavar="character"),
  make_option(c("--pct"), type="double", default=NULL, 
              help="a numeric value indicating percentage of experssion cutoff. [default= %default]", metavar="double"),
  make_option(c("--prefix"), type="character", default=NULL, 
              help="Prefix to be added to file names. [default= %default]", metavar="character")
  )


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(openxlsx)
library(dplyr)
options(stringsAsFactors = FALSE)

for(i in c("fc", "pval", "pct")) opt[[i]] <- as.numeric(opt[[i]])
geneinfo <- read.delim(opt$geneinfo, header = TRUE)
colnames(geneinfo)[c(1,2)] <- c("gene_id", "gene_name")
geneinfo <- cbind(geneinfo[,1:2], data.frame(gene_name_unique = make.unique(geneinfo$gene_name)), geneinfo[,3:ncol(geneinfo)])
message('loading tables')
de_tables <- list.files("./", pattern = paste0(opt$prefix, ".*"))
if(opt$reorder_sheets == "1") de_tables <- de_tables[order(as.numeric(gsub(".*cluster|.txt", "", basename(de_tables))),decreasing = FALSE)]
de_tables_list <- lapply(de_tables, function(i) {
  tmp <- read.delim(i, header = TRUE, row.names = NULL)
  idx_fc <- grep("log2FC$", colnames(tmp))
  if(length(idx_fc) > 1) fc_vec <- apply(tmp[,idx_fc], 1, mean) else fc_vec <- unlist(tmp[,idx_fc])
  tmp <- tmp[order(fc_vec, decreasing = TRUE),]
  colnames(tmp)[1] <- "gene_name_unique"
  tmp <- right_join(geneinfo, tmp, by = "gene_name_unique")
  })
de_tables_list_filter <- lapply(1:length(de_tables_list), function(i){
  tmp <- de_tables_list[[i]]
  if(opt$pval_flag == "1") idx_pval <- grep("p_val_adj$", colnames(tmp)) else idx_pval <- grep("p_val$", colnames(tmp))
  idx_fc <- grep("log2FC$", colnames(tmp))
  idx_pct1 <- grep("pct.1$", colnames(tmp))
  idx_pct2 <- grep("pct.2$", colnames(tmp))
  idx_de <- lapply(1:length(idx_pval), function(j) {
    if(opt$direction == "up") which(tmp[,idx_pval[j]] <= opt$pval & tmp[,idx_fc[j]] >= log2(opt$fc) & (tmp[,idx_pct1[j]] >= opt$pct | tmp[,idx_pct2[j]] >= opt$pct)) else which(tmp[,idx_pval[j]] <= opt$pval & abs(tmp[,idx_fc[j]]) >= log2(opt$fc) & (tmp[,idx_pct1[j]] >= opt$pct | tmp[,idx_pct2[j]] >= opt$pct))
    })
  idx_de <- Reduce(intersect, idx_de)
  tmp[idx_de,]
})

message(paste0('number of DEGs across clusters: ',paste0(sapply(de_tables_list_filter, function(i) nrow(i)),collapse = ",")))

message('saving tables to excel')
sheetnames <- gsub("\\.[^.]+$","",basename(de_tables))
wb <- createWorkbook()
for(i in 1:length(sheetnames)) {
  addWorksheet(wb, sheetnames[i])
  writeData(wb, sheet = sheetnames[i], de_tables_list[[i]],rowNames=FALSE,colNames=TRUE)
}
info <- data.frame(name = c("p_val_*", "avg_log2FC_*", "pct.1_*", "pct.2_*", "p_val_adj_*", "max_pval", "minimump_p_val"),
                   description = c("raw p value without adjusting for multiple hypothesis testing", "average log2FC between group1 vs group2. Positive avg_log2FC indicates upregulation in group1; negative avg_log2FC indicates downregulation in group1", "proportion of expressing cells for a certain gene in group 1", "proportion of expressing cells for a certain gene in group 2", "p value after adjusting for multiple hypothesis testing using Bonferronni method", "Largest p value of p value calculated by each group. Only applicable to FindConservedMarkers.", "Combined p value based on metap. Only applicable to FindConservedMarker."))
addWorksheet(wb, "info")
writeData(wb, sheet = "info", info, rowNames=FALSE,colNames=TRUE)
saveWorkbook(wb, paste0("./", opt$prefix,"gene_table.xlsx"), overwrite = TRUE)

wb <- createWorkbook()
for(i in 1:length(sheetnames)) {
  addWorksheet(wb, sheetnames[i])
  writeData(wb, sheet = sheetnames[i], de_tables_list_filter[[i]],rowNames=FALSE,colNames=TRUE)
}
addWorksheet(wb, "info")
writeData(wb, sheet = "info", info, rowNames=FALSE,colNames=TRUE)
saveWorkbook(wb, paste0("./", opt$prefix,"gene_table_filter.xlsx"), overwrite = TRUE)