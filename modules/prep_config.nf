process WRITECONFIGFILE {
  cpus 1
  memory '1 GB'

  publishDir(
    path: "${params.output_dir}/config/",
    mode: "copy",
    pattern: "config.txt"
  )
  
  input:
  val samplesheet
  val workflowpath
  val authorname
  val data_type
  val feature_list
  val output_dir
  val geneinfo
  val qc_only
  val ambient_RNA_removal_flag
  val doublet_removal_flag
  val adaptive_cutoff_flag
  val mt_cutoff
  val hb_cutoff
  val nFeature_cutoff
  val nCount_cutoff
  val nCell_cutoff
  val norm_dimreduc
  val norm_diff
  val cellcycle_correction_flag
  val genelist_S_phase
  val genelist_G2M_phase
  val merge_analysis
  val integration_analysis
  val merge_only
  val integration_only
  val integration_method
  val sketch_flag
  val resolution
  val vismethod
  val spatial_cluster
  val lambda
  val k_geom
  val control_var
  val case_var
  val covariate_list
  val test
  val fc
  val pval
  val pval_flag
  val pct
  
  output:
  path "config.txt"

  shell:
  """ 
  ## General
  echo "samplesheet=!{samplesheet}" > config.txt
  echo "workflowpath=!{workflowpath}" >> config.txt
  echo "authorname=!{authorname}" >> config.txt
  echo "data_type=!{data_type}" >> config.txt
  echo "feature_list=!{feature_list}" >> config.txt
  echo "output_dir=!{output_dir}" >> config.txt
  echo "geneinfo=!{geneinfo}" >> config.txt

  ## QC
  echo "qc_only=!{qc_only}" >> config.txt
  echo "ambient_RNA_removal_flag=!{ambient_RNA_removal_flag}" >> config.txt
  echo "doublet_removal_flag=!{doublet_removal_flag}" >> config.txt
  echo "adaptive_cutoff_flag=!{adaptive_cutoff_flag}" >> config.txt
  echo "mt_cutoff=!{mt_cutoff}" >> config.txt
  echo "hb_cutoff=!{hb_cutoff}" >> config.txt
  echo "nFeature_cutoff=!{nFeature_cutoff}" >> config.txt
  echo "nCount_cutoff=!{nCount_cutoff}" >> config.txt
  echo "nCell_cutoff=!{nCell_cutoff}" >> config.txt

  ## Normalization
  echo "norm_dimreduc=!{norm_dimreduc}" >> config.txt
  echo "norm_diff=!{norm_diff}" >> config.txt
  echo "cellcycle_correction_flag=!{cellcycle_correction_flag}" >> config.txt
  echo "genelist_S_phase=!{genelist_S_phase}" >> config.txt
  echo "genelist_G2M_phase=!{genelist_G2M_phase}" >> config.txt

  ## Analysis strategy
  echo "merge_analysis=!{merge_analysis}" >> config.txt
  echo "integration_analysis=!{integration_analysis}" >> config.txt
  echo "merge_only=!{merge_only}" >> config.txt
  echo "integration_only=!{integration_only}" >> config.txt

  ## Integration strategy
  echo "integration_method=!{integration_method}" >> config.txt
  echo "sketch_flag=!{sketch_flag}" >> config.txt

  ## Clustering
  echo "resolution=!{resolution}" >> config.txt
  echo "vismethod=!{vismethod}" >> config.txt
  echo "spatial_cluster=!{spatial_cluster}" >> config.txt
  echo "lambda=!{lambda}" >> config.txt
  echo "k_geom=!{k_geom}" >> config.txt 

  ## Differential expression
  echo "control_var=!{control_var}" >> config.txt
  echo "case_var=!{case_var}" >> config.txt
  echo "covariate_list=!{covariate_list}" >> config.txt
  echo "test=!{test}" >> config.txt
  echo "fc=!{fc}" >> config.txt
  echo "pval=!{pval}" >> config.txt
  echo "pval_flag=!{pval_flag}" >> config.txt
  echo "pct=!{pct}" >> config.txt
  """
}