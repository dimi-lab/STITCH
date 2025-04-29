#!/usr/bin/env nextflow


// Load modules
include { WRITECONFIGFILE } from './modules/prep_config.nf'

include { RUNQC; QCSUMMARY; LOADSAMPLE } from './modules/qc.nf'

include { SVG; DECONVOLUTION; MAPPING; INTEGRATESAMPLES; GENERATEPSEUDOBULK; FINDMARKERS; COMBINETABLES; FINALREPORT } from './modules/analysis.nf'

include { MERGESAMPLES; GENERATEPSEUDOBULK as GENERATEPSEUDOBULK_MERGE; FINDMARKERS as FINDMARKERS_MERGE; COMBINETABLES as COMBINETABLES_MERGE; FINALREPORT as FINALREPORT_MERGE } from './modules/analysis.nf'

samples_ch = Channel.fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple( row.sampleid, row.condition, file(row.secondary_output) ) }

workflow {

    WRITECONFIGFILE(params.samplesheet, params.workflowpath, params.authorname, params.data_type, params.feature_list, params.output_dir, params.geneinfo, params.parallel_strategy, params.nworkers, params.qc_only, params.binsize, params.ambient_RNA_removal_flag, params.doublet_removal_flag, params.adaptive_cutoff_flag, params.mt_cutoff, params.hb_cutoff, params.nFeature_cutoff, params.nCount_cutoff, params.nCell_cutoff, params.norm_dimreduc, params.norm_diff, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, params.svg_analysis, params.svg_method, params.deconvolution_analysis, params.reference_deconvolution, params.reference_assay_deconvolution, params.query_assay_deconvolution, params.refdata_deconvolution, params.doublet_mode, params.gene_list_reg, params.mapping_analysis, params.reference_mapping, params.reference_assay_mapping, params.query_assay_mapping, params.refdata_mapping, params.reference_reduction, params.normalization_method, params.prediction_assay, params.reduction_model, params.merge_analysis, params.integration_analysis, params.merge_only, params.integration_only, params.integration_method, params.sketch_flag, params.resolution, params.vismethod, params.spatial_cluster, params.lambda, params.k_geom, params.pseudobulk_flag, params.control_var, params.case_var, params.covariate_list, params.test, params.fc, params.pval, params.pval_flag, params.pct)
    
    RUNQC(params.workflowpath, params.data_type, params.binsize, params.ambient_RNA_removal_flag,params.doublet_removal_flag,params.adaptive_cutoff_flag,params.mt_cutoff,params.hb_cutoff,params.nFeature_cutoff,params.nCount_cutoff,params.nCell_cutoff, samples_ch)

    QCSUMMARY(params.workflowpath, params.data_type, params.ambient_RNA_removal_flag,params.doublet_removal_flag,params.adaptive_cutoff_flag,params.mt_cutoff,params.hb_cutoff,params.nFeature_cutoff,params.nCount_cutoff,params.nCell_cutoff, RUNQC.output.qc_metrics_cells.collect(), RUNQC.output.qc_metrics_summary.collect(), RUNQC.output.opt_postqc.collect(), RUNQC.output.median_umi.collect())

    if (!params.qc_only) {
        
        LOADSAMPLE(params.workflowpath, params.data_type, params.binsize, params.ambient_RNA_removal_flag,params.doublet_removal_flag,params.adaptive_cutoff_flag,params.mt_cutoff,params.hb_cutoff,params.nFeature_cutoff,params.nCount_cutoff,params.nCell_cutoff,params.geneinfo, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, QCSUMMARY.output.min_median_umi, params.norm_dimreduc, params.norm_diff,samples_ch)

        if (params.svg_analysis) {

            SVG(params.norm_dimreduc, params.svg_method, LOADSAMPLE.output.seurat_obj.collect(), samples_ch)

        }

        if (params.deconvolution_analysis) {
            
            DECONVOLUTION(params.reference_deconvolution, params.query_assay_deconvolution, params.reference_assay_deconvolution, params.refdata_deconvolution, params.doublet_mode, params.parallel_strategy, params.nworkers, params.gene_list_reg,LOADSAMPLE.output.seurat_obj.collect(), samples_ch)

        }

        if (params.mapping_analysis) {
            
            MAPPING(params.reference_mapping, params.reference_assay_mapping, params.query_assay_mapping, params.reference_reduction, params.normalization_method, params.refdata_mapping, params.prediction_assay, params.reduction_model, LOADSAMPLE.output.seurat_obj.collect(), samples_ch)

            if (params.integration_analysis) {

                INTEGRATESAMPLES(params.workflowpath, params.samplesheet,params.data_type, params.resolution,params.geneinfo, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, params.integration_method, params.sketch_flag, params.norm_dimreduc, params.spatial_cluster, params.lambda, params.k_geom, params.parallel_strategy, params.nworkers, MAPPING.output.seurat_obj.collect())

            }

            if (params.merge_analysis) {

                MERGESAMPLES(params.workflowpath, params.samplesheet,params.data_type, params.resolution,params.geneinfo, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, params.sketch_flag, params.norm_dimreduc, params.spatial_cluster, params.lambda, params.k_geom, params.parallel_strategy, params.nworkers, MAPPING.output.seurat_obj.collect())

            }

        } else {

            if (params.integration_analysis) {

                INTEGRATESAMPLES(params.workflowpath, params.samplesheet,params.data_type, params.resolution,params.geneinfo, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, params.integration_method, params.sketch_flag, params.norm_dimreduc, params.spatial_cluster, params.lambda, params.k_geom, params.parallel_strategy, params.nworkers, LOADSAMPLE.output.seurat_obj.collect())

            }

            if (params.merge_analysis) {

                MERGESAMPLES(params.workflowpath, params.samplesheet,params.data_type, params.resolution,params.geneinfo, params.cellcycle_correction_flag, params.genelist_S_phase, params.genelist_G2M_phase, params.sketch_flag, params.norm_dimreduc, params.spatial_cluster, params.lambda, params.k_geom, params.parallel_strategy, params.nworkers, LOADSAMPLE.output.seurat_obj.collect())

            }
        }

        if (params.integration_analysis && !params.integration_only) {

            cluster_ch = INTEGRATESAMPLES.output.seurat_clusters
            .splitCsv(header: true, sep: '\t')
            .map { row -> row.clusternum }
            
            on_disk_mat = INTEGRATESAMPLES.output.on_disk_mat.ifEmpty {
                empty_file = file("$projectDir/assets/NO_FILE")
                return [empty_file]
                }   

            if (params.pseudobulk_flag == "1") {

                GENERATEPSEUDOBULK(INTEGRATESAMPLES.output.integrated_obj, params.data_type,on_disk_mat)

                on_disk_mat = file("$projectDir/assets/NO_FILE")

                FINDMARKERS(params.workflowpath,params.data_type,params.pseudobulk_flag,params.control_var, params.case_var, params.covariate_list, params.test, params.sketch_flag, params.norm_diff, GENERATEPSEUDOBULK.output.seurat_obj_pb, on_disk_mat, GENERATEPSEUDOBULK.output.seurat_clusters_condition ,cluster_ch)

            } else {

            FINDMARKERS(params.workflowpath,params.data_type,params.pseudobulk_flag,params.control_var, params.case_var, params.covariate_list, params.test, params.sketch_flag, params.norm_diff, INTEGRATESAMPLES.output.integrated_obj, on_disk_mat, INTEGRATESAMPLES.output.seurat_clusters_condition ,cluster_ch)

            }

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

        if (params.merge_analysis && !params.merge_only) {

			cluster_merge_ch = MERGESAMPLES.output.seurat_clusters
            .splitCsv(header: true, sep: '\t')
            .map { row -> row.clusternum }

            on_disk_mat_merge = MERGESAMPLES.output.on_disk_mat.ifEmpty {
                empty_file = file("$projectDir/assets/NO_FILE")
                return [empty_file]
                }   

            if (params.pseudobulk_flag == "1") {

                GENERATEPSEUDOBULK_MERGE(MERGESAMPLES.output.merged_obj, params.data_type,on_disk_mat)

                on_disk_mat = file("$projectDir/assets/NO_FILE")

                FINDMARKERS_MERGE(params.workflowpath,params.data_type,params.pseudobulk_flag,params.control_var, params.case_var, params.covariate_list, params.test, params.sketch_flag, params.norm_diff, GENERATEPSEUDOBULK_MERGE.output.seurat_obj_pb, on_disk_mat, GENERATEPSEUDOBULK_MERGE.output.seurat_clusters_condition ,cluster_ch)

            } else {

            FINDMARKERS_MERGE(params.workflowpath,params.data_type,params.pseudobulk_flag,params.control_var, params.case_var, params.covariate_list, params.test, params.sketch_flag, params.norm_diff, MERGESAMPLES.output.merged_obj, on_disk_mat_merge, MERGESAMPLES.output.seurat_clusters_condition,cluster_merge_ch)
            
            }

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



