#!/usr/bin/env nextflow


// Load modules
include { WRITECONFIGFILE } from './modules/prep_config.nf'

include { RUNQC; QCSUMMARY; LOADSAMPLE } from './modules/qc.nf'

include { INTEGRATESAMPLES; FINDMARKERS; COMBINETABLES; FINALREPORT } from './modules/analysis.nf'

include { MERGESAMPLES; FINDMARKERS as FINDMARKERS_MERGE; COMBINETABLES as COMBINETABLES_MERGE; FINALREPORT as FINALREPORT_MERGE } from './modules/analysis.nf'

samples_ch = Channel.fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple( row.sampleid, row.condition, file(row.secondary_output) ) }

workflow {

    WRITECONFIGFILE(params.samplesheet, params.workflowpath, params.authorname, params.data_type, params.feature_list, params.output_dir, params.geneinfo, params.qc_only, params.ambient_RNA_removal_flag, params.doublet_removal_flag, params.adaptive_cutoff_flag, params.mt_cutoff, params.hb_cutoff, params.nFeature_cutoff, params.nCount_cutoff, params.nCell_cutoff, params.norm_dimreduc, params.norm_diff, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, params.merge_analysis, params.integration_analysis, params.merge_only, params.integration_only, params.integration_method, params.sketch_flag, params.resolution, params.vismethod, params.spatial_cluster, params.lambda, params.k_geom, params.control_var, params.case_var, params.covariate_list, params.test, params.fc, params.pval, params.pval_flag, params.pct)
    
    RUNQC(params.workflowpath, params.data_type, params.ambient_RNA_removal_flag,params.doublet_removal_flag,params.adaptive_cutoff_flag,params.mt_cutoff,params.hb_cutoff,params.nFeature_cutoff,params.nCount_cutoff,params.nCell_cutoff, samples_ch)

    QCSUMMARY(params.workflowpath, params.data_type, params.ambient_RNA_removal_flag,params.doublet_removal_flag,params.adaptive_cutoff_flag,params.mt_cutoff,params.hb_cutoff,params.nFeature_cutoff,params.nCount_cutoff,params.nCell_cutoff, RUNQC.output.qc_metrics_cells.collect(), RUNQC.output.qc_metrics_summary.collect(), RUNQC.output.opt_postqc.collect(), RUNQC.output.median_umi.collect())

    if (!params.qc_only) {
        
        LOADSAMPLE(params.workflowpath, params.data_type, params.ambient_RNA_removal_flag,params.doublet_removal_flag,params.adaptive_cutoff_flag,params.mt_cutoff,params.hb_cutoff,params.nFeature_cutoff,params.nCount_cutoff,params.nCell_cutoff,params.geneinfo, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, QCSUMMARY.output.min_median_umi, params.norm_dimreduc, params.norm_diff,samples_ch)

        if (params.integration_analysis) {

            INTEGRATESAMPLES(params.workflowpath, params.samplesheet,params.data_type, params.resolution,params.geneinfo, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, params.integration_method, params.sketch_flag, params.norm_dimreduc, params.spatial_cluster, params.lambda, params.k_geom, LOADSAMPLE.output.seurat_obj.collect())

            if (!params.integration_only) {
            cluster_ch = INTEGRATESAMPLES.output.seurat_clusters
            .splitCsv(header: true, sep: '\t')
            .map { row -> row.clusternum }

            on_disk_mat = INTEGRATESAMPLES.output.on_disk_mat.ifEmpty {
                empty_file = file("$projectDir/assets/NO_FILE")
                return [empty_file]
                }   

            FINDMARKERS(params.workflowpath,params.data_type,params.control_var, params.case_var, params.covariate_list, params.test, params.sketch_flag, params.norm_diff, INTEGRATESAMPLES.output.integrated_obj, on_disk_mat, INTEGRATESAMPLES.output.seurat_clusters_condition ,cluster_ch)

            diff_gene_tables = FINDMARKERS.output.diff_gene_table.collect().ifEmpty {
                empty_file = file("$projectDir/assets/NO_FILE")
                return [empty_file]
                }

            COMBINETABLES(params.workflowpath, params.geneinfo, params.fc, params.pval, params.pval_flag, params.pct, FINDMARKERS.output.marker_gene_table.collect(), diff_gene_tables, "integrate")

            diff_gene_table_filter = COMBINETABLES.output.diff_gene_table_filter.ifEmpty {
                empty_file = file("$projectDir/assets/NO_FILE_DIFF")
                return [empty_file]
                }

            FINALREPORT(params.workflowpath, params.data_type, params.authorname, params.samplesheet, params.sketch_flag, params.feature_list, params.cellcycle_correction_flag, params.vismethod, INTEGRATESAMPLES.output.integrated_obj, on_disk_mat, COMBINETABLES.output.marker_gene_table_filter, diff_gene_table_filter, "integrate")

            }

        }

        if (params.merge_analysis) {

            MERGESAMPLES(params.workflowpath, params.samplesheet,params.data_type, params.resolution,params.geneinfo, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, params.sketch_flag, params.norm_dimreduc, params.spatial_cluster, params.lambda, params.k_geom, LOADSAMPLE.output.seurat_obj.collect())

            if (!params.merge_only) {
			cluster_merge_ch = MERGESAMPLES.output.seurat_clusters
            .splitCsv(header: true, sep: '\t')
            .map { row -> row.clusternum }

            on_disk_mat_merge = MERGESAMPLES.output.on_disk_mat.ifEmpty {
                empty_file = file("$projectDir/assets/NO_FILE")
                return [empty_file]
                }   
                
            FINDMARKERS_MERGE(params.workflowpath,params.data_type,params.control_var, params.case_var, params.covariate_list, params.test, params.sketch_flag, params.norm_diff, MERGESAMPLES.output.merged_obj, on_disk_mat_merge, MERGESAMPLES.output.seurat_clusters_condition,cluster_merge_ch)

            diff_gene_tables_merge = FINDMARKERS_MERGE.output.diff_gene_table.collect().ifEmpty {
                empty_file = file("$projectDir/assets/NO_FILE")
                return [empty_file]
                }

            COMBINETABLES_MERGE(params.workflowpath, params.geneinfo, params.fc, params.pval, params.pval_flag, params.pct, FINDMARKERS_MERGE.output.marker_gene_table.collect(), diff_gene_tables_merge, "merge")

            diff_gene_table_filter_merge = COMBINETABLES_MERGE.output.diff_gene_table_filter.ifEmpty {
                empty_file = file("$projectDir/assets/NO_FILE_DIFF")
                return [empty_file]
                }

            FINALREPORT_MERGE(params.workflowpath, params.data_type, params.authorname, params.samplesheet, params.sketch_flag, params.feature_list, params.cellcycle_correction_flag, params.vismethod, MERGESAMPLES.output.merged_obj, on_disk_mat_merge, COMBINETABLES_MERGE.output.marker_gene_table_filter, diff_gene_table_filter_merge,"merge")
            }
        }
    }
}



