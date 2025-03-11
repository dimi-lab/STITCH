process DECONVOLUTION {
  cpus 5
  memory '50 GB'

  publishDir(
    path: "${params.output_dir}/deconvolution/",
    mode: "copy"
  )
  
  input:
  val reference_deconvolution
  val query_deconvolution
  val query_assay_deconvolution
  val reference_assay_deconvolution
  val refdata_deconvolution
  val doublet_mode
  val parallel_strategy
  val nworkers
  path seurat_obj_collect
  tuple val(sampleid), val(condition), path(secondary_output)

	output:
  path "${sampleid}_RCTD_results.rds"

	script:
	"""
  export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/deconvolution.r --reference $reference_deconvolution --query $query_deconvolution --reference_assay $reference_assay_deconvolution --query_assay $query_assay_deconvolution --refdata $refdata --doublet_mode $doublet_mode --parallel_strategy $parallel_strategy --nworkers $nworkers
	"""
}

process MAPPING {
  cpus 5
  memory '50 GB'

  input:
  val reference_mapping
  val query_mapping
  val reference_assay_mapping
  val query_assay_mapping
  val reference_reduction
  val normalization_method
  val refdata_mapping
  val prediction_assay
  val reduction_model
  path seurat_obj_collect
  tuple val(sampleid), val(condition), path(secondary_output)

  output:
  path "${sampleid}_mapped_to_ref.rds", emit: seurat_obj

  script:
  """
  export PROJECT_DIR=${projectDir}
  Rscript ${projectDir}/scripts/map_to_reference.r --reference $reference_mapping --query $query_mapping --reference_assay $reference_assay_mapping --query_assay $query_assay_mapping --reference_reduction $reference_reduction --normalization_method $normalization_method --refdata $refdata_mapping --prediction_assay $prediction_assay --reduction_model $reduction_model
  """
}

process INTEGRATESAMPLES {
  cpus 5
  memory '50 GB'

  publishDir(
    path: "${params.output_dir}/integrate/",
    mode: "copy"
  )
  
  input:
  val workflowpath
  val samplesheet
  val data_type
  val resolution
  val geneinfo
  val cellcycle_correction_flag
  val genelist_S_phase
  val genelist_G2M_phase
  val integration_method
  val sketch_flag
  val norm_dimreduc
  val spatial_cluster
  val lambda
  val k_geom
  path seurat_obj_collect

	output:
  path "integrated_obj.rds", emit: integrated_obj
	path "seurat_clusters.txt", emit: seurat_clusters
  path "seurat_clusters_condition.txt", emit: seurat_clusters_condition
  path "on_disk_mat/", emit: on_disk_mat, optional: true

	script:
	"""
  export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/integrate_samples_wrapper.r --workflowpath $workflowpath --samplesheet $samplesheet --data_type $data_type --resolution $resolution --geneinfo $geneinfo --cellcycle_correction_flag $cellcycle_correction_flag --genelist_S_phase $genelist_S_phase --genelist_G2M_phase $genelist_G2M_phase --integration_method $integration_method --sketch_flag $sketch_flag --norm_dimreduc $norm_dimreduc --spatial_cluster $spatial_cluster --lambda $lambda --k_geom $k_geom
	"""
}

process FINDMARKERS {
  cpus 1
  memory '5 GB'

  tag { "cluster_${clusternum}" }

  input:
  val workflowpath
  val data_type
  val control_var
  val case_var
  val covariate_list
  val test
  val sketch_flag
  val norm_diff
  path integrated_obj
  path on_disk_mat
  path seurat_clusters_condition
  val clusternum

	output:
  path "marker_gene_cluster${clusternum}.txt", emit: marker_gene_table
  path "diff_gene_cluster${clusternum}.txt", emit: diff_gene_table, optional: true

	script:
	"""
  export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/find_markers_wrapper.r --workflowpath $workflowpath --data_type $data_type --control_var $control_var --case_var $case_var --covariate_list $covariate_list --test $test --sketch_flag $sketch_flag --norm_diff $norm_diff --seurat_obj $integrated_obj --seurat_clusters_condition $seurat_clusters_condition --clusternum ${clusternum}
	"""
}

process COMBINETABLES {
  cpus 1
  memory '5 GB'

  publishDir(
    path: "${params.output_dir}/${integration_type}/gene_table/",
    mode: "copy",
    pattern: "*.xlsx"
  )

  input:
  val workflowpath
  val geneinfo
  val fc
  val pval
  val pval_flag
  val pct
  path marker_gene_table_collect
  path diff_gene_table_collect
  val integration_type

	output:
  path "marker_gene_table.xlsx"
  path "marker_gene_table_filter.xlsx", emit: marker_gene_table_filter
  path "diff_gene_table.xlsx", optional: true
  path "diff_gene_table_filter.xlsx", emit: diff_gene_table_filter, optional: true

	script:
	"""
  export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/combine_tables_wrapper.r --workflowpath $workflowpath --geneinfo $geneinfo --fc $fc --pval $pval --pval_flag $pval_flag --pct $pct
	"""
}

process FINALREPORT {
  cpus 1
  memory '5 GB'

  publishDir(
    path: "${params.output_dir}/${integration_type}/",
    mode: "copy",
    pattern: "Final_report.html"
  )

  publishDir(
    path: "${params.output_dir}/${integration_type}/figures/",
    mode: "copy",
    pattern: "*.png"
  )

  publishDir(
    path: "${params.output_dir}/${integration_type}/figures/",
    mode: "copy",
    pattern: "*.pdf"
  )

  publishDir(
    path: "${params.output_dir}/${integration_type}/tables/",
    mode: "copy",
    pattern: "*.txt"
  )

  input:
  val workflowpath
  val data_type
  val authorname
  val samplesheet
  val sketch_flag
  val feature_list
  val cellcycle_correction_flag
  val vismethod
  path seurat_obj
  path on_disk_mat
  path marker_gene_filtered
  path diff_gene_filtered
  val integration_type

	output:
  path "Final_report.html"
  path "*.txt"
  path "*.pdf"
  path "*.png"

	script:
	"""
  export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/final_report_wrapper.r --workflowpath $workflowpath --data_type $data_type --authorname $authorname --samplesheet $samplesheet --sketch_flag $sketch_flag --feature_list $feature_list --cellcycle_correction_flag $cellcycle_correction_flag --vismethod $vismethod --seurat_obj $seurat_obj --marker_gene_filtered $marker_gene_filtered --diff_gene_filtered $diff_gene_filtered
	"""
}

process MERGESAMPLES {
  cpus 5
  memory '50 GB'

  publishDir(
    path: "${params.output_dir}/merge/",
    mode: "copy"
  )

  input:
  val workflowpath
  val samplesheet
  val data_type
  val resolution
  val geneinfo
  val cellcycle_correction_flag
  val genelist_S_phase
  val genelist_G2M_phase
  val sketch_flag
  val norm_dimreduc
  val spatial_cluster
  val lambda
  val k_geom
  path seurat_obj_collect

	output:
  path "merged_obj.rds", emit: merged_obj
	path "seurat_clusters.txt", emit: seurat_clusters
  path "seurat_clusters_condition.txt", emit: seurat_clusters_condition
  path "on_disk_mat/", emit: on_disk_mat, optional: true

	script:
	"""
  export PROJECT_DIR=${projectDir}
	Rscript ${projectDir}/scripts/merge_samples_wrapper.r --workflowpath $workflowpath --samplesheet $samplesheet --data_type $data_type --resolution $resolution --geneinfo $geneinfo --cellcycle_correction_flag $cellcycle_correction_flag --genelist_S_phase $genelist_S_phase --genelist_G2M_phase $genelist_G2M_phase --sketch_flag $sketch_flag --norm_dimreduc $norm_dimreduc --spatial_cluster $spatial_cluster --lambda $lambda --k_geom $k_geom
	"""
}