/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPROCESS_MINFI       } from '../modules/local/dnam/preprocess_minfi/main'
include { P_VAL_CORRECTION       } from '../modules/local/dnam/p_val_correction/main'
include { FILTER_COMMON_PROBES   } from '../modules/local/dnam/filter_common_probes/main'
include { FILTER_BY_MISSING      } from '../modules/local/dnam/filter_by_missing/main'
include { FILTER_BY_VARIANCE     } from '../modules/local/dnam/filter_by_variance/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DNAM {

    take:
    ch_samplesheet
    ch_precomputed_beta_matrix
    ch_precomputed_pvals

    main:

    ch_versions = channel.empty()

    //
    // MODULE: Preprocess methylation array data with minfi
    //
    def use_precomputed = params.use_precomputed_dnam

    def ch_beta_matrix
    def ch_detection_pvals

    if (use_precomputed) {
        ch_beta_matrix = ch_precomputed_beta_matrix
        ch_detection_pvals = ch_precomputed_pvals
    } else {
        PREPROCESS_MINFI (
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(PREPROCESS_MINFI.out.versions)
        ch_beta_matrix = PREPROCESS_MINFI.out.betas
        ch_detection_pvals = PREPROCESS_MINFI.out.detection_pvals
    }

    //
    // MODULE: Replace beta values with NaN where detection p-value >= threshold
    //
    P_VAL_CORRECTION (
        ch_beta_matrix,
        ch_detection_pvals
    )
    ch_versions = ch_versions.mix(P_VAL_CORRECTION.out.versions)

    //
    // MODULE: Deduplicate probes and filter to probes common across 450K/EPIC1/EPIC2
    //
    ch_common_probes = channel.fromPath(params.common_probes, checkIfExists: true)
    FILTER_COMMON_PROBES (
        P_VAL_CORRECTION.out.corrected_betas,
        ch_common_probes
    )
    ch_versions = ch_versions.mix(FILTER_COMMON_PROBES.out.versions)

    //
    // MODULE: Remove features with more than params.missing_threshold missing values
    //
    FILTER_BY_MISSING (
        FILTER_COMMON_PROBES.out.filtered_betas
    )
    ch_versions = ch_versions.mix(FILTER_BY_MISSING.out.versions)

    //
    // MODULE: Remove features with variance <= params.variance_threshold
    //
    FILTER_BY_VARIANCE (
        FILTER_BY_MISSING.out.missing_filtered_betas
    )
    ch_versions = ch_versions.mix(FILTER_BY_VARIANCE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'methylml_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
