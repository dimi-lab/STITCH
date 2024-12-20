process RUNQC {
  cpus 1
  memory '5 GB'

  tag { "${sampleid}" }

  publishDir(
    path: "${params.output_dir}/qc/${sampleid}/figures/",
    mode: "copy",
    pattern: "*.png"
  )
  publishDir(
    path: "${params.output_dir}/qc/${sampleid}/tables/",
    mode: "copy",
    pattern: "*.txt"
  )
  publishDir(
    path: "${params.output_dir}/qc/${sampleid}/",
    mode: "copy",
    pattern: "*.html"
  )

  input:
  val workflowpath
  val data_type
  val ambient_RNA_removal_flag
  val doublet_removal_flag
  val adaptive_cutoff_flag
  val mt_cutoff
  val hb_cutoff
  val nFeature_cutoff
  val nCount_cutoff
  val nCell_cutoff
  tuple val(sampleid), val(condition), path(secondary_output)

	output:
  path "QC_report.html"
	path "qc_scatterplot.png"
  path "qc_violinplot.png"
  path "qc_summary.png"
  path "${sampleid}_qc_metrics_cells.txt", emit: qc_metrics_cells
  path "${sampleid}_qc_metrics_genes.txt"
  path "${sampleid}_qc_metrics_summary.txt", emit: qc_metrics_summary
  path "${sampleid}_opt_postQC.txt", emit: opt_postqc
  path "${sampleid}_median_umi.txt", emit: median_umi  

	script:
	"""
	export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/qc_sample_wrapper.r --workflowpath $workflowpath --data_type $data_type --ambient_RNA_removal_flag $ambient_RNA_removal_flag --doublet_removal_flag $doublet_removal_flag --adaptive_cutoff_flag $adaptive_cutoff_flag --mt_cutoff $mt_cutoff --hb_cutoff $hb_cutoff --nFeature_cutoff $nFeature_cutoff --nCount_cutoff $nCount_cutoff --nCell_cutoff $nCell_cutoff --sampleid $sampleid --condition $condition --secondary_output $secondary_output
	"""
}


process QCSUMMARY {
  cpus 1
  memory '5 GB'

  publishDir(
    path: "${params.output_dir}/qc/",
    mode: "copy",
    pattern: "*.html"
  )

  input:
  val workflowpath
  val data_type
  val ambient_RNA_removal_flag
  val doublet_removal_flag
  val adaptive_cutoff_flag
  val mt_cutoff
  val hb_cutoff
  val nFeature_cutoff
  val nCount_cutoff
  val nCell_cutoff
  path qc_metrics_cells
  path qc_metrics_summary
  path opt_postqc
  path median_umi

	output:
  path "QC_report_summary.html"
  path "min_median_umi.txt", emit: min_median_umi

	script:
	"""
	export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/qc_summary_wrapper.r --workflowpath $workflowpath --data_type $data_type --ambient_RNA_removal_flag $ambient_RNA_removal_flag --doublet_removal_flag $doublet_removal_flag --adaptive_cutoff_flag $adaptive_cutoff_flag --mt_cutoff $mt_cutoff --hb_cutoff $hb_cutoff --nFeature_cutoff $nFeature_cutoff --nCount_cutoff $nCount_cutoff --nCell_cutoff $nCell_cutoff
	"""
}

process LOADSAMPLE {
  cpus 1
  memory '5 GB'

  tag { "${sampleid}" }

  input:
  val workflowpath
  val data_type
  val ambient_RNA_removal_flag
  val doublet_removal_flag
  val adaptive_cutoff_flag
  val mt_cutoff
  val hb_cutoff
  val nFeature_cutoff
  val nCount_cutoff
  val nCell_cutoff
  val geneinfo
  val cellcycle_correction_flag
  val genelist_S_phase
  val genelist_G2M_phase
  path min_median_umi
  val norm_dimreduc
  val norm_diff
  tuple val(sampleid), val(condition), path(secondary_output)

	output:
  path "${sampleid}.rds", emit: seurat_obj

	script:
	"""
	export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/data_loader_wrapper.r --workflowpath $workflowpath --data_type $data_type --ambient_RNA_removal_flag $ambient_RNA_removal_flag --doublet_removal_flag $doublet_removal_flag --adaptive_cutoff_flag $adaptive_cutoff_flag --mt_cutoff $mt_cutoff --hb_cutoff $hb_cutoff --nFeature_cutoff $nFeature_cutoff --nCount_cutoff $nCount_cutoff --nCell_cutoff $nCell_cutoff --geneinfo $geneinfo --cellcycle_correction_flag $cellcycle_correction_flag --genelist_S_phase $genelist_S_phase --genelist_G2M_phase $genelist_G2M_phase --min_median_umi $min_median_umi --norm_dimreduc $norm_dimreduc --norm_diff $norm_diff --sampleid $sampleid --condition $condition --secondary_output $secondary_output 
	"""
}