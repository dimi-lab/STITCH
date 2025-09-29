loweroutlier_IQR <- function(vec1) {
  tmp <- quantile(vec1,na.rm=T)
  as.numeric(tmp[2] - 1.5*(tmp[4] - tmp[2]))}

higheroutlier_IQR <- function(vec1) {
  tmp <- quantile(vec1,na.rm=T)
  as.numeric(tmp[4] + 1.5*(tmp[4] - tmp[2]))}

remove_object <- function(object_name){
  for(i in object_name) assign(i, NULL,envir = .GlobalEnv)
  rm(list = object_name, envir = .GlobalEnv)
  invisible(gc(verbose = FALSE))
}

#' @param sampleid Character. Sample ID. Default is NULL.
#' @param condition Character. Condition. Default is NULL.
#' @param secondary_output Character. Path to secondary output. Default is NULL.
#' @param stage Character. Either "QC" or "loader", indicating it is QC or loader step. Default is "QC".
#' @param workflowpath Character. Full path to workflow directory. Default is NULL.
#' @param binsize bin size to read in. Default is NULL.
#' @param adaptive_cutoff_flag Character 0 or 1 to indicate whether to apply adaptive cutoff identification (based on IQR). Default is NULL.
#' @param nCount_cutoff Double. Cutoff for the total number of UMI counts. Default is 200.
#' @param nFeature_cutoff Double. Cutoff for the total number of detectable genes/features. Default is 50.
#' @param mt_cutoff Double. Cutoff for the percentage of mitochondria concentration, e.g., 40 indicates 40%. Default is 40.
#' @param hb_cutoff Double. Cutoff for the percentage of hemoglobin concentration, e.g., 20 indicates 20%. Default is 20.
#' @param nCell_cutoff Double. Cutoff for the number of cells with expression for a feature/gene. Default is 10.
#' @param geneinfo Character. Path to gene-level annotation file. This is used to add feature-level metadata. Only applicable for the 'loader' stage. Default is "NA".
#' @param cellcycle_correction_flag Character 0 or 1 indicating whether to estimate and correct for cell-cycle effects. If set to 1, also need to specify gene lists for S and G2M phases. Default is "1".
#' @param genelist_S_phase Character. Path to gene list (Gene Symbols) for cell-cycle S-phase, one gene per line. Default is "NA".
#' @param genelist_G2M_phase Character. Path to gene list (Gene Symbols) for cell-cycle G2M-phase, one gene per line. Default is "NA".
#' @param min_median_umi Character. Path to min_median_umi file. Default is "NA".
#' @param norm_dimreduc Character. normalization method for dimension reduction. Default is "NA".
#' @param norm_diff Character. normalization method for differential testin. Default is "NA".
#' @return seurat object

qcsample_VisiumHD <- function(
    sampleid = NULL,
    condition = NULL,
    secondary_output = NULL,
    stage = "QC",
    workflowpath = NULL,
    binsize = NULL,
    adaptive_cutoff_flag = NULL,
    nCount_cutoff = 200,
    nFeature_cutoff = 50,
    mt_cutoff = 40,
    hb_cutoff = 20,
    nCell_cutoff = 10,
    geneinfo = "NA",
    cellcycle_correction_flag = "1",
    genelist_S_phase = "NA",
    genelist_G2M_phase = "NA",
    min_median_umi = "NA",
    norm_dimreduc = "NA",
    norm_diff = "NA"
) {
  sampleinfo <- data.frame(sampleid = sampleid, condition = condition, secondary_output = secondary_output)
  
  seurat_obj <- Load10X_Spatial(data.dir = sampleinfo$secondary_output, slice = sampleid, bin.size = as.integer(binsize))
  metadata <- sampleinfo[rep(1,ncol(seurat_obj)),]
  rownames(metadata) <- colnames(seurat_obj)
  seurat_obj <- AddMetaData(seurat_obj, metadata)
  
  message(paste0("calculate mt percent for sample "),sampleid)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-|^mt-")
  message(paste0("calculate HB percent for sample "),sampleid)
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^HB.*|^hb.*")
  idx_nCount_Spatial <- grep("nCount_Spatial", colnames(seurat_obj@meta.data))
  idx_nFeature_Spatial <- grep("nFeature_Spatial", colnames(seurat_obj@meta.data))
  assay <- Assays(seurat_obj)
  if(stage == "QC") {
    write.table(data.frame(sampleid = seurat_obj$sampleid,
                           barcode = colnames(seurat_obj),
                           nCount_Spatial = seurat_obj@meta.data[,idx_nCount_Spatial],
                           nFeature_Spatial = seurat_obj@meta.data[,idx_nFeature_Spatial],
                           percent.mt = seurat_obj$percent.mt,
                           percent.hb = seurat_obj$percent.hb),paste0(sampleid,'_qc_metrics_cells.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=F)}
  if(adaptive_cutoff_flag == "1"){
    nCount_cutoff <- max(nCount_cutoff,loweroutlier_IQR(seurat_obj@meta.data[,idx_nCount_Spatial]))
    nFeature_cutoff <- max(nFeature_cutoff,loweroutlier_IQR(seurat_obj@meta.data[,idx_nFeature_Spatial]))
    mt_cutoff <- min(mt_cutoff,higheroutlier_IQR(seurat_obj$percent.mt))
    hb_cutoff <- min(hb_cutoff,higheroutlier_IQR(seurat_obj$percent.hb))
  }
  idx <- which(seurat_obj@meta.data[,idx_nCount_Spatial] >= nCount_cutoff & seurat_obj@meta.data[,idx_nFeature_Spatial] >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff & seurat_obj$percent.hb <= hb_cutoff)
  gene_summary <- data.frame(GeneName = rownames(seurat_obj),
                             nCell = apply(GetAssayData(seurat_obj, assay = assay, layer = "counts")[,idx],1,function(i) sum(i>0)))
  if(stage == "QC") write.table(gene_summary,paste0(sampleid,'_qc_metrics_genes.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=F)
  if(adaptive_cutoff_flag == "1") nCell_cutoff <- max(nCell_cutoff,loweroutlier_IQR(gene_summary$nCell))
  idx1 <- which(gene_summary$nCell >= nCell_cutoff)
  if(stage == "QC"){
    median_umi <- median(apply(GetAssayData(seurat_obj, assay = assay, layer = "counts")[idx1,idx],2,sum))
    write.table(median_umi,paste0(sampleid,'_median_umi.txt'),sep='\t',row.names=FALSE,col.names=FALSE,quote=F)
    qc_summary <- data.frame(sampleid = sampleid,
                             total_bins = ncol(seurat_obj),
                             nCount_cutoff = sum(seurat_obj@meta.data[,idx_nCount_Spatial] >= nCount_cutoff),
                             nFeature_cutoff = sum(seurat_obj@meta.data[,idx_nCount_Spatial] >= nCount_cutoff & seurat_obj@meta.data[,idx_nFeature_Spatial] >= nFeature_cutoff),
                             mt_cutoff = sum(seurat_obj@meta.data[,idx_nCount_Spatial] >= nCount_cutoff & seurat_obj@meta.data[,idx_nFeature_Spatial] >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff),
                             hb_cutoff = sum(seurat_obj@meta.data[,idx_nCount_Spatial] >= nCount_cutoff & seurat_obj@meta.data[,idx_nFeature_Spatial] >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff & seurat_obj$percent.hb <= hb_cutoff))
    qc_summary$total_genes <- sum(gene_summary$nCell >= nCell_cutoff)
    write.table(qc_summary,paste0(sampleid,'_qc_metrics_summary.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(list(sampleid = sampleid, 
                     nCount_cutoff = nCount_cutoff,
                     nFeature_cutoff = nFeature_cutoff,
                     mt_cutoff = mt_cutoff,
                     hb_cutoff = hb_cutoff,
                     nCell_cutoff = nCell_cutoff,
                     adaptive_cutoff_flag = adaptive_cutoff_flag), paste0(sampleid,"_opt_postQC.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
  } else if(stage == "loader"){
    seurat_obj <- subset(seurat_obj, cells = idx, features = idx1)
    message(paste0('predict cell cycle phase for sample ', sampleid))
    ### use normalize data before cellcyclescoring
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, assay = assay)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
    seurat_obj <- ScaleData(seurat_obj)
    if(geneinfo != "NA"){
      gene_metadata <- data.table::fread(geneinfo, header = TRUE,data.table = FALSE)
      colnames(gene_metadata)[1:2] <- c("gene_id", "gene_name")
      gene_metadata$gene_name_unique <- make.unique(gene_metadata$gene_name)
      rownames(gene_metadata) <- gene_metadata$gene_name_unique
      seurat_obj[[assay]] <- AddMetaData(seurat_obj[[assay]], metadata = gene_metadata)
    }
    if(cellcycle_correction_flag == "1") {
      s.features <- read.delim(genelist_S_phase, header = FALSE)$V1
      g2m.features <- read.delim(genelist_G2M_phase, header = FALSE)$V1
      seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.features, g2m.features = g2m.features, set.ident = FALSE)
      vars.to.regress <- c('percent.mt','Phase')
    } else  vars.to.regress <- c('percent.mt')
    if(norm_dimreduc == "SCT" | norm_diff == "SCT"){
      message(paste0('perform SCTransform for sample ',sampleid))
      min_median_umi <- read.delim(min_median_umi, header = FALSE)$V1
      seurat_obj <- SCTransform(seurat_obj, method="glmGamPoi", assay = assay, variable.features.n = 3000,vst.flavor = "v2",vars.to.regress = vars.to.regress,verbose = FALSE,return.only.var.genes = FALSE, scale_factor=min_median_umi)
      if(geneinfo != "NA") seurat_obj@assays$SCT <- AddMetaData(seurat_obj@assays$SCT, metadata = gene_metadata)
    }
    seurat_obj <- AddMetaData(seurat_obj, metadata = GetTissueCoordinates(seurat_obj))
  }
  return(seurat_obj)
}

#' @param sampleid Character. Sample ID. Default is NULL.
#' @param condition Character. Condition. Default is NULL.
#' @param secondary_output Character. Path to secondary output. Default is NULL.
#' @param stage Character. Either "QC" or "loader", indicating it is QC or loader step. Default is "QC".
#' @param workflowpath Character. Full path to workflow directory. Default is NULL.
#' @param adaptive_cutoff_flag Character 0 or 1 to indicate whether to apply adaptive cutoff identification (based on IQR). Default is NULL.
#' @param nCount_cutoff Double. Cutoff for the total number of UMI counts. Default is 200.
#' @param nFeature_cutoff Double. Cutoff for the total number of detectable genes/features. Default is 50.
#' @param mt_cutoff Double. Cutoff for the percentage of mitochondria concentration, e.g., 40 indicates 40%. Default is 40.
#' @param hb_cutoff Double. Cutoff for the percentage of hemoglobin concentration, e.g., 20 indicates 20%. Default is 20.
#' @param nCell_cutoff Double. Cutoff for the number of cells with expression for a feature/gene. Default is 10.
#' @param geneinfo Character. Path to gene-level annotation file. This is used to add feature-level metadata. Only applicable for the 'loader' stage. Default is "NA".
#' @param cellcycle_correction_flag Character 0 or 1 indicating whether to estimate and correct for cell-cycle effects. If set to 1, also need to specify gene lists for S and G2M phases. Default is "1".
#' @param genelist_S_phase Character. Path to gene list (Gene Symbols) for cell-cycle S-phase, one gene per line. Default is "NA".
#' @param genelist_G2M_phase Character. Path to gene list (Gene Symbols) for cell-cycle G2M-phase, one gene per line. Default is "NA".
#' @param min_median_umi Character. Path to min_median_umi file. Default is "NA".
#' @param norm_dimreduc Character. normalization method for dimension reduction. Default is "NA".
#' @param norm_diff Character. normalization method for differential testin. Default is "NA".
#' @return seurat object

qcsample_Visium <- function(
    sampleid = NULL,
    condition = NULL,
    secondary_output = NULL,
    stage = "QC",
    workflowpath = NULL,
    adaptive_cutoff_flag = NULL,
    nCount_cutoff = 200,
    nFeature_cutoff = 50,
    mt_cutoff = 40,
    hb_cutoff = 20,
    nCell_cutoff = 10,
    geneinfo = "NA",
    cellcycle_correction_flag = "1",
    genelist_S_phase = "NA",
    genelist_G2M_phase = "NA",
    min_median_umi = "NA",
    norm_dimreduc = "NA",
    norm_diff = "NA"
) {
  sampleinfo <- data.frame(sampleid = sampleid, condition = condition, secondary_output = secondary_output)
  
  seurat_obj <- Load10X_Spatial(data.dir = sampleinfo$secondary_output, slice = sampleid)
  metadata <- sampleinfo[rep(1,ncol(seurat_obj)),]
  rownames(metadata) <- colnames(seurat_obj)
  seurat_obj <- AddMetaData(seurat_obj, metadata)
  
  message(paste0("calculate mt percent for sample "),sampleid)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-|^mt-")
  message(paste0("calculate HB percent for sample "),sampleid)
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^HB.*|^hb.*")
  if(stage == "QC") {
    write.table(data.frame(sampleid = seurat_obj$sampleid,
                           barcode = colnames(seurat_obj),
                           nCount_Spatial = seurat_obj$nCount_Spatial,
                           nFeature_Spatial = seurat_obj$nFeature_Spatial,
                           percent.mt = seurat_obj$percent.mt,
                           percent.hb = seurat_obj$percent.hb),paste0(sampleid,'_qc_metrics_cells.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=F)}
  if(adaptive_cutoff_flag == "1"){
    nCount_cutoff <- max(nCount_cutoff,loweroutlier_IQR(seurat_obj$nCount_Spatial))
    nFeature_cutoff <- max(nFeature_cutoff,loweroutlier_IQR(seurat_obj$nFeature_Spatial))
    mt_cutoff <- min(mt_cutoff,higheroutlier_IQR(seurat_obj$percent.mt))
    hb_cutoff <- min(hb_cutoff,higheroutlier_IQR(seurat_obj$percent.hb))
  }
  idx <- which(seurat_obj$nCount_Spatial >= nCount_cutoff & seurat_obj$nFeature_Spatial >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff & seurat_obj$percent.hb <= hb_cutoff)
  gene_summary <- data.frame(GeneName = rownames(seurat_obj),
                             nCell = apply(GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")[,idx],1,function(i) sum(i>0)))
  if(stage == "QC") write.table(gene_summary,paste0(sampleid,'_qc_metrics_genes.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=F)
  if(adaptive_cutoff_flag == "1") nCell_cutoff <- max(nCell_cutoff,loweroutlier_IQR(gene_summary$nCell))
  idx1 <- which(gene_summary$nCell >= nCell_cutoff)
  if(stage == "QC"){
    median_umi <- median(apply(GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")[idx1,idx],2,sum))
    write.table(median_umi,paste0(sampleid,'_median_umi.txt'),sep='\t',row.names=FALSE,col.names=FALSE,quote=F)
    qc_summary <- data.frame(sampleid = sampleid,
                             total_spots = ncol(seurat_obj),
                             nCount_cutoff = sum(seurat_obj$nCount_Spatial >= nCount_cutoff),
                             nFeature_cutoff = sum(seurat_obj$nCount_Spatial >= nCount_cutoff & seurat_obj$nFeature_Spatial >= nFeature_cutoff),
                             mt_cutoff = sum(seurat_obj$nCount_Spatial >= nCount_cutoff & seurat_obj$nFeature_Spatial >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff),
                             hb_cutoff = sum(seurat_obj$nCount_Spatial >= nCount_cutoff & seurat_obj$nFeature_Spatial >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff & seurat_obj$percent.hb <= hb_cutoff))
    qc_summary$total_genes <- sum(gene_summary$nCell >= nCell_cutoff)
    write.table(qc_summary,paste0(sampleid,'_qc_metrics_summary.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(list(sampleid = sampleid, 
                     nCount_cutoff = nCount_cutoff,
                     nFeature_cutoff = nFeature_cutoff,
                     mt_cutoff = mt_cutoff,
                     hb_cutoff = hb_cutoff,
                     nCell_cutoff = nCell_cutoff,
                     adaptive_cutoff_flag = adaptive_cutoff_flag), paste0(sampleid,"_opt_postQC.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
  } else if(stage == "loader"){
    seurat_obj <- subset(seurat_obj, cells = idx, features = idx1)
    message(paste0('predict cell cycle phase for sample ', sampleid))
    ### use normalize data before cellcyclescoring
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, assay = "Spatial")
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
    seurat_obj <- ScaleData(seurat_obj)
    if(geneinfo != "NA"){
      gene_metadata <- data.table::fread(geneinfo, header = TRUE,data.table = FALSE)
      colnames(gene_metadata)[1:2] <- c("gene_id", "gene_name")
      gene_metadata$gene_name_unique <- make.unique(gene_metadata$gene_name)
      rownames(gene_metadata) <- gene_metadata$gene_name_unique
      seurat_obj@assays$Spatial <- AddMetaData(seurat_obj@assays$Spatial, metadata = gene_metadata)
    }
    if(cellcycle_correction_flag == "1") {
      s.features <- read.delim(genelist_S_phase, header = FALSE)$V1
      g2m.features <- read.delim(genelist_G2M_phase, header = FALSE)$V1
      seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.features, g2m.features = g2m.features, set.ident = FALSE)
      vars.to.regress <- c('percent.mt','Phase')
    } else  vars.to.regress <- c('percent.mt')
    if(norm_dimreduc == "SCT" | norm_diff == "SCT"){
      message(paste0('perform SCTransform for sample ',sampleid))
      min_median_umi <- read.delim(min_median_umi, header = FALSE)$V1
      seurat_obj <- SCTransform(seurat_obj, method="glmGamPoi", assay = "Spatial", variable.features.n = 3000,vst.flavor = "v2",vars.to.regress = vars.to.regress,verbose = FALSE,return.only.var.genes = FALSE, scale_factor=min_median_umi)
      if(geneinfo != "NA") seurat_obj@assays$SCT <- AddMetaData(seurat_obj@assays$SCT, metadata = gene_metadata)
    }
    seurat_obj <- AddMetaData(seurat_obj, metadata = GetTissueCoordinates(seurat_obj))
  }
  return(seurat_obj)
}


#' @param sampleid Character. Sample ID. Default is NULL.
#' @param condition Character. Condition. Default is NULL.
#' @param secondary_output Character. Path to secondary output. Default is NULL.
#' @param stage Character. Either "QC" or "loader", indicating it is QC or loader step. Default is "QC".
#' @param workflowpath Character. Full path to workflow directory. Default is NULL.
#' @param ambient_RNA_removal_flag Character. 0 or 1 to indicate whether to perform ambient RNA removal/correction using SoupX. Default is NULL.
#' @param doublet_removal_flag Character. 0 or 1 to indicate whether to perform doublet removal using scDblFinder. Default is NULL.
#' @param adaptive_cutoff_flag Character 0 or 1 to indicate whether to apply adaptive cutoff identification (based on IQR). Default is NULL.
#' @param nCount_cutoff Double. Cutoff for the total number of UMI counts. Default is 200.
#' @param nFeature_cutoff Double. Cutoff for the total number of detectable genes/features. Default is 50.
#' @param mt_cutoff Double. Cutoff for the percentage of mitochondria concentration, e.g., 40 indicates 40%. Default is 40.
#' @param hb_cutoff Double. Cutoff for the percentage of hemoglobin concentration, e.g., 20 indicates 20%. Default is 20.
#' @param nCell_cutoff Double. Cutoff for the number of cells with expression for a feature/gene. Default is 10.
#' @param geneinfo Character. Path to gene-level annotation file. This is used to add feature-level metadata. Only applicable for the 'loader' stage. Default is "NA".
#' @param cellcycle_correction_flag Character. 0 or 1 indicating whether to estimate and correct for cell-cycle effects. If set to 1, also need to specify gene lists for S and G2M phases. Default is "1".
#' @param genelist_S_phase Character. Path to gene list (Gene Symbols) for cell-cycle S-phase, one gene per line. Default is "NA".
#' @param genelist_G2M_phase Character. Path to gene list (Gene Symbols) for cell-cycle G2M-phase, one gene per line. Default is "NA".
#' @param min_median_umi Character. Path to min_median_umi file. Default is "NA".
#' @param norm_dimreduc Character. normalization method for dimension reduction. Default is "NA".
#' @param norm_diff Character. normalization method for differential testin. Default is "NA".
#' @return seurat object

qcsample_scRNAseq <- function(
    sampleid = NULL,
    condition = NULL,
    secondary_output = NULL,
    stage = "QC",
    workflowpath = NULL,
    ambient_RNA_removal_flag = NULL,
    doublet_removal_flag = NULL,
    adaptive_cutoff_flag = NULL,
    nCount_cutoff = 200,
    nFeature_cutoff = 50,
    mt_cutoff = 40,
    hb_cutoff = 20,
    nCell_cutoff = 10,
    geneinfo = "NA",
    cellcycle_correction_flag = "1",
    genelist_S_phase = "NA",
    genelist_G2M_phase = "NA",
    min_median_umi = "NA",
    norm_dimreduc = "NA",
    norm_diff = "NA"
) {
  sampleinfo <- data.frame(sampleid = sampleid, condition = condition, secondary_output = secondary_output)
  if(length(list.files(sampleinfo[[3]], pattern = "*raw_feature_bc_matrix$"))) raw_feature_bc_matrix_dir <- file.path(sampleinfo[[3]], list.files(sampleinfo[[3]], pattern = "*raw_feature_bc_matrix$")) else raw_feature_bc_matrix_dir <- "NA"
  
  message(paste0('create seurat object for ', sampleid))
  data <- Read10X(data.dir=file.path(sampleinfo[[3]], list.files(sampleinfo[[3]], pattern = "*filtered_feature_bc_matrix$")), gene.column = 2)
  message(paste0("read10X completed for sample ",sampleid))
  metainfo <- sapply(sampleinfo,function(i) rep(i,ncol(data)))
  rownames(metainfo) <- colnames(data)
  colnames(metainfo)[1:3] <- c('sampleid','condition','secondary_output')
  seurat_obj <- CreateSeuratObject(counts=data, project=sampleid,meta.data = as.data.frame(metainfo),min.cells=0,min.features=0)
  remove_object("data")
  message("perform first round data normalization using LogNormalize")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(object = seurat_obj, assay = "RNA", verbose = FALSE ,npcs = 30)
  seurat_obj <- FindNeighbors(object = seurat_obj, assay = "RNA", reduction = "pca", dims = 1:30, verbose = FALSE)
  ## todo: change to Leiden algorithm
  message("perform unsupervised clustering using Louvain algorithm")
  seurat_obj <- FindClusters(object = seurat_obj, algorithm = 1, verbose = FALSE, cluster.name = "seurat_clusters")
  clusters <- seurat_obj$seurat_clusters
  meta.data <- seurat_obj@meta.data
  toc <- GetAssayData(seurat_obj, assay = "RNA",layer = "count")
  remove_object("seurat_obj")
  if(stage == "QC" | (stage == "loader" & ambient_RNA_removal_flag == "1")){
    message("Perform ambient RNA estimation using SoupX")
    if(raw_feature_bc_matrix_dir != "NA"){
      message("raw_feature_bc_matrix exists, use droplets to profile Soup")
      tod <- Read10X(raw_feature_bc_matrix_dir, gene.column = 2)
      tod <- tod[rownames(tod) %in% rownames(toc),]
      sc <- SoupChannel(tod, toc)
    } else {
      message("raw_feature_bc_matrix not present, use filtered droplets to estimate Soup")
      sc <- SoupChannel(toc, toc, calcSoupProfile = FALSE)
      soupProf <- data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
      sc <- setSoupProfile(sc, soupProf)
    }
    remove_object("tod")
    sc <- setClusters(sc, clusters = clusters)
    tmp <- try(autoEstCont(sc,doPlot = FALSE), silent = TRUE)
    if(class(tmp)[1] == "try-error") {
      message("auto estimation failed for SoupX, set contamination fraction to 0.1")
      sc <- setContaminationFraction(sc, 0.1)
    } else sc <- autoEstCont(sc,doPlot = FALSE)
    if(ambient_RNA_removal_flag == "0") out <- toc else out <- adjustCounts(sc, roundToInt = TRUE, clusters = clusters)
    rho <- sc$metaData$rho
    remove_object(c("sc", "toc"))
    seurat_obj <- CreateSeuratObject(counts=out, project=sampleid,meta.data = as.data.frame(meta.data),min.cells=0,min.features=0)
    seurat_obj[["soupX.rho"]] <- rho
    remove_object(c("out"))
  } else {
    seurat_obj <- CreateSeuratObject(counts=toc, project=sampleid,meta.data = as.data.frame(meta.data),min.cells=0,min.features=0)
    remove_object("toc")
  }
  remove_object("meta.data")
  if(stage == "QC" | (stage == "loader" & doublet_removal_flag == "1")){
    message("perform doublet detection using scDblFinder")
    if(length(unique(clusters)) ==1) sce <- scDblFinder(GetAssayData(seurat_obj, assay = "RNA",layer = "count"), cluster = 2) else sce <- scDblFinder(GetAssayData(seurat_obj, assay = "RNA",layer = "count"), clusters=clusters)
    seurat_obj[["scDblFinder.score"]] <- sce$scDblFinder.score
    seurat_obj[["scDblFinder.class"]] <- sce$scDblFinder.class
  }
  message(paste0("calculate mt percent for sample "),sampleid)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-|^mt-")
  message(paste0("calculate HB percent for sample "),sampleid)
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^HB.*|^hb.*")
  if(stage == "QC") {
    write.table(data.frame(sampleid = seurat_obj$sampleid,
                           barcode = colnames(seurat_obj),
                           nCount_RNA = seurat_obj$nCount_RNA,
                           nFeature_RNA = seurat_obj$nFeature_RNA,
                           percent.mt = seurat_obj$percent.mt,
                           percent.hb = seurat_obj$percent.hb,
                           scDblFinder.score = seurat_obj$scDblFinder.score,
                           scDblFinder.class = seurat_obj$scDblFinder.class,
                           soupX.rho = seurat_obj$soupX.rho
                           ),paste0(sampleid,'_qc_metrics_cells.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=F)}
  if(adaptive_cutoff_flag == "1"){
    nCount_cutoff <- max(nCount_cutoff,loweroutlier_IQR(seurat_obj$nCount_RNA))
    nFeature_cutoff <- max(nFeature_cutoff,loweroutlier_IQR(seurat_obj$nFeature_RNA))
    mt_cutoff <- min(mt_cutoff,higheroutlier_IQR(seurat_obj$percent.mt))
    hb_cutoff <- min(hb_cutoff,higheroutlier_IQR(seurat_obj$percent.hb))
  }
  if(doublet_removal_flag == "1") idx <- which(seurat_obj$scDblFinder.class != "doublet" & seurat_obj$nCount_RNA >= nCount_cutoff & seurat_obj$nFeature_RNA >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff & seurat_obj$percent.hb <= hb_cutoff) else idx <- which(seurat_obj$nCount_RNA >= nCount_cutoff & seurat_obj$nFeature_RNA >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff & seurat_obj$percent.hb <= hb_cutoff)
  
  gene_summary <- data.frame(GeneName = rownames(seurat_obj),
                             nCell = apply(GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[,idx],1,function(i) sum(i>0)))
  if(stage == "QC") write.table(gene_summary,paste0(sampleid,'_qc_metrics_genes.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=F)
  if(adaptive_cutoff_flag == "1") nCell_cutoff <- max(nCell_cutoff,loweroutlier_IQR(gene_summary$nCell))
  idx1 <- which(gene_summary$nCell >= nCell_cutoff)
  if(stage == "QC"){
    median_umi <- median(apply(GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[idx1,idx],2,sum))
    write.table(median_umi,paste0(sampleid,'_median_umi.txt'),sep='\t',row.names=FALSE,col.names=FALSE,quote=F)
    if(doublet_removal_flag == "1"){
      qc_summary <- data.frame(sampleid = sampleid,
                               total_cells = ncol(seurat_obj),
                               doublet_removal_flag = sum(seurat_obj$scDblFinder.class != "doublet"),
                               nCount_cutoff = sum(seurat_obj$scDblFinder.class != "doublet" & seurat_obj$nCount_RNA >= nCount_cutoff),
                               nFeature_cutoff = sum(seurat_obj$scDblFinder.class != "doublet" & seurat_obj$nCount_RNA >= nCount_cutoff & seurat_obj$nFeature_RNA >= nFeature_cutoff),
                               mt_cutoff = sum(seurat_obj$scDblFinder.class != "doublet" & seurat_obj$nCount_RNA >= nCount_cutoff & seurat_obj$nFeature_RNA >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff),
                               hb_cutoff = sum(seurat_obj$scDblFinder.class != "doublet" & seurat_obj$nCount_RNA >= nCount_cutoff & seurat_obj$nFeature_RNA >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff & seurat_obj$percent.hb <= hb_cutoff))
    } else {
      qc_summary <- data.frame(sampleid = sampleid,
                               total_cells = ncol(seurat_obj),
                               nCount_cutoff = sum(seurat_obj$nCount_RNA >= nCount_cutoff),
                               nFeature_cutoff = sum(seurat_obj$nCount_RNA >= nCount_cutoff & seurat_obj$nFeature_RNA >= nFeature_cutoff),
                               mt_cutoff = sum(seurat_obj$nCount_RNA >= nCount_cutoff & seurat_obj$nFeature_RNA >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff),
                               hb_cutoff = sum(seurat_obj$nCount_RNA >= nCount_cutoff & seurat_obj$nFeature_RNA >= nFeature_cutoff & seurat_obj$percent.mt <= mt_cutoff & seurat_obj$percent.hb <= hb_cutoff))
    }
    qc_summary$total_genes <- sum(gene_summary$nCell >= nCell_cutoff)
    write.table(qc_summary,paste0(sampleid,'_qc_metrics_summary.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(list(sampleid = sampleid,
                     ambient_RNA_removal_flag = ambient_RNA_removal_flag,
                     doublet_removal_flag = doublet_removal_flag,
                     nCount_cutoff = nCount_cutoff,
                     nFeature_cutoff = nFeature_cutoff,
                     mt_cutoff = mt_cutoff,
                     hb_cutoff = hb_cutoff,
                     nCell_cutoff = nCell_cutoff,
                     adaptive_cutoff_flag = adaptive_cutoff_flag), paste0(sampleid,"_opt_postQC.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
    } else if(stage == "loader"){
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[idx], features = rownames(seurat_obj)[idx1])
    message(paste0('predict cell cycle phase for sample ', sampleid))
    ### use normalize data before cellcyclescoring
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)
    if(geneinfo != "NA"){
      gene_metadata <- data.table::fread(geneinfo, header = TRUE,data.table = FALSE)
      colnames(gene_metadata)[1:2] <- c("gene_id", "gene_name")
      gene_metadata$gene_name_unique <- make.unique(gene_metadata$gene_name)
      rownames(gene_metadata) <- gene_metadata$gene_name_unique
      seurat_obj@assays$RNA <- AddMetaData(seurat_obj@assays$RNA, metadata = gene_metadata)
    }
    if(cellcycle_correction_flag == "1") {
      s.features <- read.delim(genelist_S_phase, header = FALSE)$V1
      g2m.features <- read.delim(genelist_G2M_phase, header = FALSE)$V1
      seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.features, g2m.features = g2m.features, set.ident = FALSE)
      vars.to.regress <- c('percent.mt','Phase')
    } else  vars.to.regress <- c('percent.mt')
    if(norm_dimreduc == "SCT" | norm_diff == "SCT"){
      message(paste0('perform SCTransform for sample ',sampleid))
      min_median_umi <- read.delim(min_median_umi, header = FALSE)$V1
      seurat_obj <- SCTransform(seurat_obj, method="glmGamPoi", assay = "RNA", variable.features.n = 3000,vst.flavor = "v2",vars.to.regress = vars.to.regress,verbose = FALSE,return.only.var.genes = FALSE, scale_factor=min_median_umi)
      if(geneinfo != "NA") seurat_obj@assays$SCT <- AddMetaData(seurat_obj@assays$SCT, metadata = gene_metadata)
    }}
  return(seurat_obj)
}












#### fix rare scenarios where NAs were generated combined p values
#### not used for now
FindConservedMarkers_fix <- function (object, ident.1, ident.2 = NULL, grouping.var, assay = "RNA", 
                                      slot = "data", meta.method = metap::minimump, verbose = TRUE, 
                                      ...) 
{
  metap.installed <- PackageCheck("metap", error = FALSE)
  if (!metap.installed[1]) {
    stop("Please install the metap package to use FindConservedMarkers.", 
         "\nThis can be accomplished with the following commands: ", 
         "\n----------------------------------------", "\ninstall.packages('BiocManager')", 
         "\nBiocManager::install('multtest')", "\ninstall.packages('metap')", 
         "\n----------------------------------------", call. = FALSE)
  }
  if (!is.function(x = meta.method)) {
    stop("meta.method should be a function from the metap package. Please see https://cran.r-project.org/web/packages/metap/metap.pdf for a detailed description of the available functions.")
  }
  object.var <- FetchData(object = object, vars = grouping.var)
  object <- SetIdent(object = object, cells = colnames(x = object), 
                     value = paste(Idents(object = object), object.var[, 
                                                                       1], sep = "_"))
  levels.split <- names(x = sort(x = table(object.var[, 1])))
  num.groups <- length(levels.split)
  cells <- list()
  for (i in 1:num.groups) {
    cells[[i]] <- rownames(x = object.var[object.var[, 1] == 
                                            levels.split[i], , drop = FALSE])
  }
  marker.test <- list()
  ident.2.save <- ident.2
  for (i in 1:num.groups) {
    level.use <- levels.split[i]
    ident.use.1 <- paste(ident.1, level.use, sep = "_")
    ident.use.1.exists <- ident.use.1 %in% Idents(object = object)
    if (!all(ident.use.1.exists)) {
      bad.ids <- ident.1[!ident.use.1.exists]
      warning("Identity: ", paste(bad.ids, collapse = ", "), 
              " not present in group ", level.use, ". Skipping ", 
              level.use, call. = FALSE, immediate. = TRUE)
      next
    }
    ident.2 <- ident.2.save
    cells.1 <- WhichCells(object = object, idents = ident.use.1)
    if (is.null(x = ident.2)) {
      cells.2 <- setdiff(x = cells[[i]], y = cells.1)
      ident.use.2 <- names(x = which(x = table(Idents(object = object)[cells.2]) > 
                                       0))
      ident.2 <- gsub(pattern = paste0("_", level.use), 
                      replacement = "", x = ident.use.2)
      if (length(x = ident.use.2) == 0) {
        stop(paste("Only one identity class present:", 
                   ident.1))
      }
    }
    else {
      ident.use.2 <- paste(ident.2, level.use, sep = "_")
    }
    if (verbose) {
      message("Testing group ", level.use, ": (", paste(ident.1, 
                                                        collapse = ", "), ") vs (", paste(ident.2, collapse = ", "), 
              ")")
    }
    ident.use.2.exists <- ident.use.2 %in% Idents(object = object)
    if (!all(ident.use.2.exists)) {
      bad.ids <- ident.2[!ident.use.2.exists]
      warning("Identity: ", paste(bad.ids, collapse = ", "), 
              " not present in group ", level.use, ". Skipping ", 
              level.use, call. = FALSE, immediate. = TRUE)
      next
    }
    marker.test[[i]] <- FindMarkers(object = object, assay = assay, 
                                    slot = slot, ident.1 = ident.use.1, ident.2 = ident.use.2, 
                                    verbose = verbose, ...)
    names(x = marker.test)[i] <- levels.split[i]
  }
  marker.test <- Filter(f = Negate(f = is.null), x = marker.test)
  genes.conserved <- Reduce(f = intersect, x = lapply(X = marker.test, 
                                                      FUN = function(x) {
                                                        return(rownames(x = x))
                                                      }))
  markers.conserved <- list()
  for (i in 1:length(x = marker.test)) {
    markers.conserved[[i]] <- marker.test[[i]][genes.conserved, 
    ]
    colnames(x = markers.conserved[[i]]) <- paste(names(x = marker.test)[i], 
                                                  colnames(x = markers.conserved[[i]]), sep = "_")
  }
  markers.combined <- Reduce(cbind, markers.conserved)
  pval.codes <- colnames(x = markers.combined)[grepl(pattern = "*_p_val$", 
                                                     x = colnames(x = markers.combined))]
  if (length(x = pval.codes) > 1) {
    markers.combined$max_pval <- apply(X = markers.combined[, 
                                                            pval.codes, drop = FALSE], MARGIN = 1, FUN = max)
    ### handle NAs in pvals
    tmp <- markers.combined[,pval.codes, drop = FALSE]
    idx <- which(apply(tmp,1,function(i) any(is.na(i))))
    tmp[is.na(tmp)] <- 1
    combined.pval <- data.frame(cp = apply(X = tmp, MARGIN = 1, FUN = function(x) {
                                                                  return(meta.method(x)$p)
                                                                }))
    combined.pval$cp[idx] <- NA
    meta.method.name <- as.character(x = formals()$meta.method)
    if (length(x = meta.method.name) == 3) {
      meta.method.name <- meta.method.name[3]
    }
    colnames(x = combined.pval) <- paste0(meta.method.name, 
                                          "_p_val")
    markers.combined <- cbind(markers.combined, combined.pval)
    markers.combined <- markers.combined[order(markers.combined[, 
                                                                paste0(meta.method.name, "_p_val")]), ]
  }
  else {
    warning("Only a single group was tested", call. = FALSE, 
            immediate. = TRUE)
  }
  return(markers.combined)
}


PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package.installed], collapse = ', '),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}


InsertColumn <- function(df1,df2,after){
  for(i in 1:length(after)){
    idx <- which(colnames(df1) == after[i])
    df1 <- cbind(df1[,1:idx],df2[i],df1[,(idx+1):ncol(df1)])
  }
  return(df1)
}
