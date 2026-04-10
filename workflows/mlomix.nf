/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GEX_REF_PREPROCESSOR } from './gex_ref_preprocessor'
include { GEX } from './gex'
include { DNAM } from './dnam'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MLOMIX {

    take:
    ch_gex_samplesheet
    ch_datasets
    ch_batches
    ch_classes
    ch_dnam_samplesheet
    ch_dnam_beta_matrix
    ch_dnam_pvals
    genome
    annotation_version
    random_seed
    run_gex
    run_dnam

    main:

    ch_versions = channel.empty()
    ch_gex_transposed = channel.empty()
    ch_classes_tsv = channel.empty()

    if (run_gex) {
        GEX_REF_PREPROCESSOR (
            genome,
            annotation_version
        )

        GEX (
            ch_gex_samplesheet,
            ch_datasets,
            ch_batches,
            ch_classes,
            GEX_REF_PREPROCESSOR.out.filtered_annotations,
            random_seed
        )

        ch_versions = ch_versions.mix(GEX_REF_PREPROCESSOR.out.versions)
        ch_versions = ch_versions.mix(GEX.out.versions)
        ch_gex_transposed = GEX.out.transposed_csv
        ch_classes_tsv = GEX.out.classes_tsv
    }

    if (run_dnam) {
        DNAM (
            ch_dnam_samplesheet,
            ch_dnam_beta_matrix,
            ch_dnam_pvals
        )
        ch_versions = ch_versions.mix(DNAM.out.versions)
    }

    emit:
    versions = ch_versions
    gex_transposed_csv = ch_gex_transposed
    classes_tsv = ch_classes_tsv

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
