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
    use_precomputed_dnam

    main:

    ch_versions = channel.empty()
    ch_gex_transposed = channel.empty()
    ch_classes_tsv = channel.empty()

    ch_genome_for_gex = run_gex
        .filter { flag -> flag }
        .map { flag -> genome }

    ch_annotation_for_gex = run_gex
        .filter { flag -> flag }
        .map { flag -> annotation_version }

    ch_gex_samplesheet_gated = ch_gex_samplesheet
        .combine(run_gex)
        .filter { item -> item[-1] }
        .map { item -> [item[0], item[1]] }

    ch_datasets_gated = ch_datasets
        .combine(run_gex)
        .filter { item -> item[-1] }
        .map { item -> [item[0], item[1], item[2]] }

    ch_batches_gated = ch_batches
        .combine(run_gex)
        .filter { item -> item[-1] }
        .map { item -> item[0] }

    ch_classes_gated = ch_classes
        .combine(run_gex)
        .filter { item -> item[-1] }
        .map { item -> item[0] }

    ch_random_seed_gated = run_gex
        .filter { flag -> flag }
        .map { flag -> random_seed }

    GEX_REF_PREPROCESSOR (
        ch_genome_for_gex,
        ch_annotation_for_gex
    )

    GEX (
        ch_gex_samplesheet_gated,
        ch_datasets_gated,
        ch_batches_gated,
        ch_classes_gated,
        GEX_REF_PREPROCESSOR.out.filtered_annotations,
        ch_random_seed_gated
    )

    ch_versions = ch_versions.mix(GEX_REF_PREPROCESSOR.out.versions)
    ch_versions = ch_versions.mix(GEX.out.versions)
    ch_gex_transposed = GEX.out.transposed_csv
    ch_classes_tsv = GEX.out.classes_tsv

    ch_dnam_samplesheet_gated = ch_dnam_samplesheet
        .combine(run_dnam)
        .filter { item -> item[-1] }
        .map { item -> item[0] }

    ch_use_precomputed_dnam_gated = use_precomputed_dnam
        .combine(run_dnam)
        .filter { item -> item[-1] }
        .map { item -> item[0] }

    DNAM (
        ch_dnam_samplesheet_gated,
        ch_dnam_beta_matrix,
        ch_dnam_pvals,
        ch_use_precomputed_dnam_gated
    )
    ch_versions = ch_versions.mix(DNAM.out.versions)

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
