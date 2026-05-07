/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { CONCATENATE_GEX        } from '../modules/local/gex/concatenate_gex'
include { UMAP as UMAP_GEX_RAW_BY_BATCH  } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_BC_BY_BATCH   } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_NORM_BY_BATCH } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_RAW_BY_CLASS  } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_BC_BY_CLASS   } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_NORM_BY_CLASS } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_PROCESSED_BY_BATCH } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_PROCESSED_BY_CLASS } from '../modules/local/umap/umap'
include { BATCH_CORRECT          } from '../modules/local/gex/batch_correct'
include { FILTER_BY_VARIANCE     } from '../modules/local/filter_by_variance/main'
include { FILTER_GENES           } from '../modules/local/gex/filter_genes/filter_genes'
include { MERGE_DATASETS         } from '../modules/local/gex/merge_datasets'
include { NORMALIZE              } from '../modules/local/gex/normalize'
include { TRANSPOSE              } from '../modules/local/transpose'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GEX {
    take:
    _ch_samplesheet       // channel: samplesheet from the samplesheet read in from --input
    ch_dataset            // channel: datasets with sample names and paths
    ch_batches            // channel: batch file information
    ch_classes            // channel: classes from samplesheet metadata
    ch_annotations        // channel: gene annotations from the GEX_REF_PREPROCESSOR subworkflow
    random_seed           //  int: Random seed for reproducibility

    main:

    ch_versions = channel.empty()
    ch_visuals = channel.empty()

    CONCATENATE_GEX (
        ch_dataset
    )
    ch_versions = ch_versions.mix(CONCATENATE_GEX.out.versions)

    FILTER_GENES (
        CONCATENATE_GEX.out.dataset_name,
        CONCATENATE_GEX.out.concatenated_gex_csv,
        params.gex_genes_file,
        ch_annotations,
    )
    ch_versions = ch_versions.mix(FILTER_GENES.out.versions)

    MERGE_DATASETS (
        FILTER_GENES.out.filtered_genes_csv.collect()
    )
    ch_versions = ch_versions.mix(MERGE_DATASETS.out.versions)

    if (params.gex_intermediate_umaps) {
        UMAP_GEX_RAW_BY_BATCH (
            "gex.raw.by_batch",
            MERGE_DATASETS.out.merged_csv,
            ch_batches,
            false,
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_GEX_RAW_BY_BATCH.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_GEX_RAW_BY_BATCH.out.umap_svg)

        UMAP_GEX_RAW_BY_CLASS (
            "gex.raw.by_class",
            MERGE_DATASETS.out.merged_csv,
            ch_classes,
            false,
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_GEX_RAW_BY_CLASS.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_GEX_RAW_BY_CLASS.out.umap_svg)
    }

    BATCH_CORRECT (
        MERGE_DATASETS.out.merged_csv,
        ch_batches,
        ch_classes
    )
    ch_versions = ch_versions.mix(BATCH_CORRECT.out.versions)

    if (params.gex_intermediate_umaps) {
        UMAP_GEX_BC_BY_BATCH (
            "gex.batch_corrected.by_batch",
            BATCH_CORRECT.out.batch_corrected_csv,
            ch_batches,
            false,
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_GEX_BC_BY_BATCH.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_GEX_BC_BY_BATCH.out.umap_svg)

        UMAP_GEX_BC_BY_CLASS (
            "gex.batch_corrected.by_class",
            BATCH_CORRECT.out.batch_corrected_csv,
            ch_classes,
            false,
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_GEX_BC_BY_CLASS.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_GEX_BC_BY_CLASS.out.umap_svg)
    }

    NORMALIZE (
        BATCH_CORRECT.out.batch_corrected_csv,
        ch_annotations
    )
    ch_versions = ch_versions.mix(NORMALIZE.out.versions)

    if (params.gex_intermediate_umaps) {
        UMAP_GEX_NORM_BY_BATCH (
            "gex.normalized.by_batch",
            NORMALIZE.out.normalized_csv,
            ch_batches,
            false,
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_GEX_NORM_BY_BATCH.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_GEX_NORM_BY_BATCH.out.umap_svg)

        UMAP_GEX_NORM_BY_CLASS (
            "gex.normalized.by_class",
            NORMALIZE.out.normalized_csv,
            ch_classes,
            false,
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_GEX_NORM_BY_CLASS.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_GEX_NORM_BY_CLASS.out.umap_svg)
    }

    // Variance filtering on normalized GEX before transpose
    FILTER_BY_VARIANCE (
        NORMALIZE.out.normalized_csv.map { f -> [ 'merged_datasets', 'merged_datasets', f ] },
        'gex'
    )
    ch_versions = ch_versions.mix(FILTER_BY_VARIANCE.out.versions)
    ch_visuals = ch_visuals.mix(FILTER_BY_VARIANCE.out.variance_plot_before_png)
    ch_visuals = ch_visuals.mix(FILTER_BY_VARIANCE.out.variance_plot_before_svg)
    ch_visuals = ch_visuals.mix(FILTER_BY_VARIANCE.out.variance_plot_after_png)
    ch_visuals = ch_visuals.mix(FILTER_BY_VARIANCE.out.variance_plot_after_svg)

    UMAP_GEX_PROCESSED_BY_BATCH (
        "gex.processed.by_batch",
        FILTER_BY_VARIANCE.out.variance_filtered_betas.map { _dataset_name, _sample_name, f -> f },
        ch_batches,
        false,
        random_seed
    )
    ch_versions = ch_versions.mix(UMAP_GEX_PROCESSED_BY_BATCH.out.versions)
    ch_visuals = ch_visuals.mix(UMAP_GEX_PROCESSED_BY_BATCH.out.umap_svg)

    UMAP_GEX_PROCESSED_BY_CLASS (
        "gex.processed.by_class",
        FILTER_BY_VARIANCE.out.variance_filtered_betas.map { _dataset_name, _sample_name, f -> f },
        ch_classes,
        false,
        random_seed
    )
    ch_versions = ch_versions.mix(UMAP_GEX_PROCESSED_BY_CLASS.out.versions)
    ch_visuals = ch_visuals.mix(UMAP_GEX_PROCESSED_BY_CLASS.out.umap_svg)

    TRANSPOSE (
        FILTER_BY_VARIANCE.out.variance_filtered_betas.map { _dataset_name, _sample_name, f -> f },
        "gex"
    )
    ch_versions = ch_versions.mix(TRANSPOSE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'gexml_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    transposed_csv = TRANSPOSE.out.transposed_csv
    classes_tsv    = ch_classes                  // channel: path(classes.tsv)
    visuals        = ch_visuals                  // channel: path(*.png|*.svg)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
