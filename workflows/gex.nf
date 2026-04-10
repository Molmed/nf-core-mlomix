/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { CLASS_REPORTER         } from '../modules/local/gex/class_reporter/class_reporter'
include { CONCATENATE_GEX        } from '../modules/local/gex/concatenate_gex'
include { UMAP as UMAP_RAW_BY_BATCH  } from '../modules/local/gex/umap/umap'
include { UMAP as UMAP_BC_BY_BATCH   } from '../modules/local/gex/umap/umap'
include { UMAP as UMAP_NORM_BY_BATCH } from '../modules/local/gex/umap/umap'
include { UMAP as UMAP_RAW_BY_CLASS  } from '../modules/local/gex/umap/umap'
include { UMAP as UMAP_BC_BY_CLASS   } from '../modules/local/gex/umap/umap'
include { UMAP as UMAP_NORM_BY_CLASS } from '../modules/local/gex/umap/umap'
include { BATCH_CORRECT          } from '../modules/local/gex/batch_correct'
include { FILTER_GENES           } from '../modules/local/gex/filter_genes/filter_genes'
include { MERGE_DATASETS         } from '../modules/local/gex/merge_datasets'
include { NORMALIZE              } from '../modules/local/gex/normalize'
include { TRANSPOSE              } from '../modules/local/gex/transpose'

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

    CLASS_REPORTER (
        ch_classes
    )
    ch_versions = ch_versions.mix(CLASS_REPORTER.out.versions)

    CONCATENATE_GEX (
        ch_dataset
    )
    ch_versions = ch_versions.mix(CONCATENATE_GEX.out.versions)

    FILTER_GENES (
        CONCATENATE_GEX.out.dataset_name,
        CONCATENATE_GEX.out.concatenated_gex_csv,
        ch_annotations,
    )
    ch_versions = ch_versions.mix(FILTER_GENES.out.versions)

    MERGE_DATASETS (
        FILTER_GENES.out.filtered_genes_csv.collect()
    )
    ch_versions = ch_versions.mix(MERGE_DATASETS.out.versions)

    UMAP_RAW_BY_BATCH (
        "raw.by_batch",
        MERGE_DATASETS.out.merged_csv,
        ch_batches,
        false,
        random_seed
    )
    ch_versions = ch_versions.mix(UMAP_RAW_BY_BATCH.out.versions)

    UMAP_RAW_BY_CLASS (
        "raw.by_class",
        MERGE_DATASETS.out.merged_csv,
        ch_classes,
        false,
        random_seed
    )
    ch_versions = ch_versions.mix(UMAP_RAW_BY_CLASS.out.versions)

    BATCH_CORRECT (
        MERGE_DATASETS.out.merged_csv,
        ch_batches,
        ch_classes
    )
    ch_versions = ch_versions.mix(BATCH_CORRECT.out.versions)

    UMAP_BC_BY_BATCH (
        "batch_corrected.by_batch",
        BATCH_CORRECT.out.batch_corrected_csv,
        ch_batches,
        false,
        random_seed
    )
    ch_versions = ch_versions.mix(UMAP_BC_BY_BATCH.out.versions)

    UMAP_BC_BY_CLASS (
        "batch_corrected.by_class",
        BATCH_CORRECT.out.batch_corrected_csv,
        ch_classes,
        false,
        random_seed
    )
    ch_versions = ch_versions.mix(UMAP_BC_BY_CLASS.out.versions)

    NORMALIZE (
        BATCH_CORRECT.out.batch_corrected_csv,
        ch_annotations
    )
    ch_versions = ch_versions.mix(NORMALIZE.out.versions)

    UMAP_NORM_BY_BATCH (
        "normalized.by_batch",
        NORMALIZE.out.normalized_csv,
        ch_batches,
        false,
        random_seed
    )
    ch_versions = ch_versions.mix(UMAP_NORM_BY_BATCH.out.versions)

    UMAP_NORM_BY_CLASS (
        "normalized.by_class",
        NORMALIZE.out.normalized_csv,
        ch_classes,
        false,
        random_seed
    )
    ch_versions = ch_versions.mix(UMAP_NORM_BY_CLASS.out.versions)

    TRANSPOSE (
        NORMALIZE.out.normalized_csv
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
