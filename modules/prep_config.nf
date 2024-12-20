process WRITECONFIGFILE {
  cpus 1
  memory '1 GB'

  publishDir(
    path: "${params.output_dir}/config/",
    mode: "copy",
    pattern: "*.txt"
  )
  
  input:
  val workflowpath
  val authorname
  val samplesheet
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
  val cellcycle_correction_flag
  val genelist_S_phase
  val genelist_G2M_phase
  val merge_analysis
  val integration_analysis
  val integration_method
  val sketch_flag
  val resolution
  val vismethod
  val control_var
  val case_var
  val covariate_list
  val test
  val fc
  val pval
  val pval_flag
  val pct
  
  output:
  path "config_QC.txt", emit: configfile_QC
  path "config_loader.txt", emit: configfile_loader
  path "config_clustering.txt", emit: configfile_clustering
  path "config_markers.txt", emit: configfile_markers
  path "config_tables.txt", emit: configfile_tables
  path "config_report.txt", emit: configfile_report

  shell:
  """ 
  ## QC config
  echo "workflowpath=!{workflowpath}" > config_QC.txt
  echo "data_type=!{data_type}" >> config_QC.txt
  echo "ambient_RNA_removal_flag=!{ambient_RNA_removal_flag}" >> config_QC.txt
  echo "doublet_removal_flag=!{doublet_removal_flag}" >> config_QC.txt
  echo "adaptive_cutoff_flag=!{adaptive_cutoff_flag}" >> config_QC.txt
  echo "mt_cutoff=!{mt_cutoff}" >> config_QC.txt
  echo "hb_cutoff=!{hb_cutoff}" >> config_QC.txt
  echo "nFeature_cutoff=!{nFeature_cutoff}" >> config_QC.txt
  echo "nCount_cutoff=!{nCount_cutoff}" >> config_QC.txt
  echo "nCell_cutoff=!{nCell_cutoff}" >> config_QC.txt

  ## loader config
  cp config_QC.txt config_loader.txt
  echo "geneinfo=!{geneinfo}" >> config_loader.txt
  echo "cellcycle_correction_flag=!{cellcycle_correction_flag}" >> config_loader.txt
  echo "genelist_S_phase=!{genelist_S_phase}" >> config_loader.txt
  echo "genelist_G2M_phase=!{genelist_G2M_phase}" >> config_loader.txt

  ## clustering config
  echo "workflowpath=!{workflowpath}" > config_clustering.txt
  echo "samplesheet=!{samplesheet}" >> config_clustering.txt
  echo "data_type=!{data_type}" >> config_clustering.txt
  echo "resolution=!{resolution}" >> config_clustering.txt
  echo "geneinfo=!{geneinfo}" >> config_clustering.txt
  echo "cellcycle_correction_flag=!{cellcycle_correction_flag}" >> config_clustering.txt
  echo "genelist_S_phase=!{genelist_S_phase}" >> config_clustering.txt
  echo "genelist_G2M_phase=!{genelist_G2M_phase}" >> config_clustering.txt
  echo "integration_method=!{integration_method}" >> config_clustering.txt
  echo "sketch_flag=!{sketch_flag}" >> config_clustering.txt

  ## findmarkers config
  echo "workflowpath=!{workflowpath}" > config_markers.txt
  echo "data_type=!{data_type}" >> config_markers.txt
  echo "control_var=!{control_var}" >> config_markers.txt
  echo "case_var=!{case_var}" >> config_markers.txt
  echo "covariate_list=!{covariate_list}" >> config_markers.txt
  echo "test=!{test}" >> config_markers.txt

  ## combine_table config
  echo "workflowpath=!{workflowpath}" > config_tables.txt
  echo "geneinfo=!{geneinfo}" >> config_tables.txt
  echo "fc=!{fc}" >> config_tables.txt
  echo "pval=!{pval}" >> config_tables.txt
  echo "pval_flag=!{pval_flag}" >> config_tables.txt
  echo "pct=!{pct}" >> config_tables.txt

  ## report config
  echo "workflowpath=!{workflowpath}" > config_report.txt
  echo "data_type=!{data_type}" >> config_report.txt
  echo "authorname=!{authorname}" >> config_report.txt
  echo "samplesheet=!{samplesheet}" >> config_report.txt
  echo "sketch_flag=!{sketch_flag}" >> config_report.txt
  echo "feature_list=!{feature_list}" >> config_report.txt
  echo "cellcycle_correction_flag=!{cellcycle_correction_flag}" >> config_report.txt
  echo "vismethod=!{vismethod}" >> config_report.txt
  """
}