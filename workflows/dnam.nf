/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPROCESS_MINFI       } from '../modules/local/dnam/preprocess_minfi/main'
include { CONCATENATE_DNAM       } from '../modules/local/dnam/concatenate_dnam'
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
    _ch_precomputed_beta_matrix
    _ch_precomputed_pvals
    _ch_use_precomputed_dnam

    main:

    ch_versions = channel.empty()

    //
    // MODULE: Concatenate precomputed DNAm matrices by sample
    //
    ch_precomputed_rows = ch_samplesheet
        .flatMap { rows -> (rows instanceof List) ? rows : [rows] }
        .filter { row -> row['dnam_beta_matrix_file'] || row['dnam_pvals_file'] }

    ch_precomputed_rows
        .map { row -> [row['dataset'], row] }
        .groupTuple()
        .map { dataset ->
            [
                dataset[0],
                dataset[1].collect { sample -> sample['id'] },
                dataset[1].collect { sample -> file(sample['dnam_beta_matrix_file'], checkIfExists: true).toString() },
                dataset[1].collect { sample -> sample['dnam_pvals_file'] ? file(sample['dnam_pvals_file'], checkIfExists: true).toString() : null }
            ]
        }
        .set { ch_precomputed_dnam }

    CONCATENATE_DNAM (
        ch_precomputed_dnam
    )
    ch_versions = ch_versions.mix(CONCATENATE_DNAM.out.versions)

    //
    // MODULE: Preprocess methylation array data with minfi
    //
    ch_samplesheet_for_minfi = ch_samplesheet
        .flatMap { rows -> (rows instanceof List) ? rows : [rows] }
        .filter { row -> row['sentrix_id'] || row['sentrix_position'] || row['idats_basename'] }
        .map { row -> "${row['id']},${row['sentrix_id']},${row['sentrix_position']},${row['idats_basename']}" }
        .collect()
        .map { lines -> "sample,sentrix_id,sentrix_position,idats_basename\n${lines.join('\n')}\n" }
        .collectFile(
            storeDir: "${params.outdir}/dnam",
            name: 'idat_samplesheet.csv',
            newLine: false
        )

    PREPROCESS_MINFI (
        ch_samplesheet_for_minfi
    )
    ch_versions = ch_versions.mix(PREPROCESS_MINFI.out.versions)

    ch_precomputed_beta_by_dataset = CONCATENATE_DNAM.out.beta_matrix
        .map { beta ->
            def dataset_name = beta.baseName.replaceFirst(/\.beta_matrix$/, '')
            [dataset_name, beta]
        }

    ch_precomputed_pvals_by_dataset = CONCATENATE_DNAM.out.detection_pvals
        .map { detection_pvals ->
            def dataset_name = detection_pvals.baseName.replaceFirst(/\.detection_pvals$/, '')
            [dataset_name, detection_pvals]
        }

    ch_precomputed_pairs = ch_precomputed_beta_by_dataset
        .join(ch_precomputed_pvals_by_dataset)

    ch_precomputed_beta_only = ch_precomputed_beta_by_dataset
        .join(ch_precomputed_pvals_by_dataset, remainder: true)
        .filter { tuple_item -> tuple_item[2] == null }
        .map { tuple_item -> tuple_item[1] }

    ch_beta_for_correction = PREPROCESS_MINFI.out.betas.mix(
        ch_precomputed_pairs.map { tuple_item -> tuple_item[1] }
    )
    ch_detection_for_correction = PREPROCESS_MINFI.out.detection_pvals.mix(
        ch_precomputed_pairs.map { tuple_item -> tuple_item[2] }
    )

    //
    // MODULE: Replace beta values with NaN where detection p-value >= threshold
    //
    P_VAL_CORRECTION (
        ch_beta_for_correction,
        ch_detection_for_correction
    )
    ch_versions = ch_versions.mix(P_VAL_CORRECTION.out.versions)

    ch_corrected_or_passthrough_betas = P_VAL_CORRECTION.out.corrected_betas.mix(ch_precomputed_beta_only)

    //
    // MODULE: Deduplicate probes and filter to probes common across 450K/EPIC1/EPIC2
    //
    ch_common_probes = channel.fromPath(params.common_probes, checkIfExists: true)
    FILTER_COMMON_PROBES (
        ch_corrected_or_passthrough_betas,
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
