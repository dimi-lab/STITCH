## Generate annotation file using biomart
## use R4.4.1 version to avoid issue with biomaRt
library(biomaRt)

## human 2024A
human_genes <- as.data.frame(rtracklayer::import("refdata-gex-GRCh38-2024-A/genes/genes.gtf", format = "gff"))
human_genes_sel <- human_genes[human_genes$type == "gene",c("gene_id","gene_name","gene_type","gene_version","seqnames", "start", "end", "width", "strand")]
listEnsemblArchives()
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "https://jul2023.archive.ensembl.org")
results <- getBM(attributes = c('ensembl_gene_id',"hgnc_symbol", "description"),
                 filters = "ensembl_gene_id", values = human_genes_sel$gene_id,
                 mart = ensembl)
human_genes_sel$gene_description <- results$description[match(human_genes_sel$gene_id,results$ensembl_gene_id)]
human_genes_sel <- human_genes_sel[,c("gene_id","gene_name","gene_type","gene_description","gene_version","seqnames", "start", "end", "width", "strand")]
write.table(human_genes_sel, "../docs/human_gene_info_2024A.tsv", row.names = F, col.names = T, sep="\t", quote = F)

## human 2020A
human_genes <- as.data.frame(rtracklayer::import("refdata-gex-GRCh38-2020-A/genes/genes.gtf", format = "gff"))
human_genes_sel <- human_genes[human_genes$type == "gene",c("gene_id","gene_name","gene_type","gene_version","seqnames", "start", "end", "width", "strand")]
listEnsemblArchives()
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "https://sep2019.archive.ensembl.org")
results <- getBM(attributes = c('ensembl_gene_id',"hgnc_symbol", "description"),
                 filters = "ensembl_gene_id", values = human_genes_sel$gene_id,
                 mart = ensembl)
human_genes_sel$gene_description <- results$description[match(human_genes_sel$gene_id,results$ensembl_gene_id)]
human_genes_sel <- human_genes_sel[,c("gene_id","gene_name","gene_type","gene_description","gene_version","seqnames", "start", "end", "width", "strand")]
write.table(human_genes_sel, "../docs/human_gene_info_2020A.tsv", row.names = F, col.names = T, sep="\t", quote = F)

## mouse 2024A
mouse_genes <- as.data.frame(rtracklayer::import("refdata-gex-GRCm39-2024-A/genes/genes.gtf", format = "gff"))
mouse_genes_sel <- mouse_genes[mouse_genes$type == "gene",c("gene_id","gene_name","gene_type","gene_version","seqnames", "start", "end", "width", "strand")]
ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl", host = "https://jul2023.archive.ensembl.org")
results <- getBM(attributes = c('ensembl_gene_id',"mgi_symbol", "mgi_description"),
                 filters = "ensembl_gene_id", values = mouse_genes_sel$gene_id,
                 mart = ensembl)
## remove duplicate records for ensembl_id
idx <- c(which(results$ensembl_gene_id == "ENSMUSG00000115016" & results$mgi_symbol == "AA536875")
)
results <- results[-idx,]
mouse_genes_sel$gene_description <- results$mgi_description[match(mouse_genes_sel$gene_id,results$ensembl_gene_id)]
mouse_genes_sel <- mouse_genes_sel[,c("gene_id","gene_name","gene_type","gene_description","gene_version","seqnames", "start", "end", "width", "strand")]
write.table(mouse_genes_sel, "../docs/mouse_gene_info_2024A.tsv", row.names = F, col.names = T, sep="\t", quote = F)

## mouse 2020A
mouse_genes <- as.data.frame(rtracklayer::import("refdata-gex-mm10-2020-A/genes/genes.gtf", format = "gff"))
mouse_genes_sel <- mouse_genes[mouse_genes$type == "gene",c("gene_id","gene_name","gene_type","gene_version","seqnames", "start", "end", "width", "strand")]
ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl", host = "https://sep2019.archive.ensembl.org")
results <- getBM(attributes = c('ensembl_gene_id',"mgi_symbol", "mgi_description"),
                 filters = "ensembl_gene_id", values = mouse_genes_sel$gene_id,
                 mart = ensembl)
## remove duplicate records for ensembl_id
idx <- c(which(results$ensembl_gene_id == "ENSMUSG00000115016" & results$mgi_symbol == "AA536875"),
         which(results$ensembl_gene_id == "ENSMUSG00000094121" & results$mgi_symbol == "Ccl21c"),
         which(results$ensembl_gene_id == "ENSMUSG00000079774" & results$mgi_symbol %in% c("Fam205a2", "Fam205a3")),
         which(results$ensembl_gene_id == "ENSMUSG00000096271" & results$mgi_symbol == "Ccl21b"),
         which(results$ensembl_gene_id == "ENSMUSG00000096873" & results$mgi_symbol == "Ccl21b")
)
results <- results[-idx,]
mouse_genes_sel$gene_description <- results$mgi_description[match(mouse_genes_sel$gene_id,results$ensembl_gene_id)]
mouse_genes_sel <- mouse_genes_sel[,c("gene_id","gene_name","gene_type","gene_description","gene_version","seqnames", "start", "end", "width", "strand")]
write.table(mouse_genes_sel, "../docs/mouse_gene_info_2020A.tsv", row.names = F, col.names = T, sep="\t", quote = F)

#### cell-cycle gene lists
human_gene_info <- read.delim("../docs/human_gene_info_2024A.tsv", header = TRUE)
write.table(cc.genes$s.genes, "../docs/S_genes_human.tsv", row.names = F, col.names = F, sep = "\t", quote = FALSE)
write.table(cc.genes$g2m.genes, "../docs/G2M_genes_human.tsv", row.names = F, col.names = F, sep = "\t", quote = FALSE)

library("biomaRt")
## use archive version to avoid issue with getLDS
human <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org", mirror = "useast")
mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org", mirror = "useast")
#### retrive information between two linked dataset
s.genes.mouse <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = cc.genes$s.genes , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
g2m.genes.mouse <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = cc.genes$g2m.genes , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
#mouse_gene_info <- read.delim("../docs/mouse_gene_info_2024A.tsv", header = TRUE)
write.table(s.genes.mouse$MGI.symbol, "../docs/S_genes_mouse.tsv", row.names = F, col.names = F, sep = "\t", quote = FALSE)
write.table(g2m.genes.mouse$MGI.symbol, "../docs/G2M_genes_mouse.tsv", row.names = F, col.names = F, sep = "\t", quote = FALSE)


