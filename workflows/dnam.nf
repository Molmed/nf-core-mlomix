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
    // Build per-sample DNAm channels from precomputed inputs
    //
    ch_precomputed_rows = ch_samplesheet
        .flatMap { rows -> (rows instanceof List) ? rows : [rows] }
        .filter { row -> row['dnam_beta_matrix_file'] }

    ch_precomputed_pairs = ch_precomputed_rows
        .filter { row -> row['dnam_pvals_file'] }
        .map { row ->
            def dataset_name = row['dataset'] == null ? 'dataset' : row['dataset'].toString().trim()
            dataset_name = dataset_name ? dataset_name : 'dataset'
            dataset_name = dataset_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            def sample_name = row['id'] == null ? 'sample' : row['id'].toString().trim()
            sample_name = sample_name ? sample_name : 'sample'
            sample_name = sample_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            [
                dataset_name,
                sample_name,
                file(row['dnam_beta_matrix_file'], checkIfExists: true),
                file(row['dnam_pvals_file'], checkIfExists: true)
            ]
        }

    ch_precomputed_beta_only = ch_precomputed_rows
        .filter { row -> !row['dnam_pvals_file'] }
        .map { row ->
            def dataset_name = row['dataset'] == null ? 'dataset' : row['dataset'].toString().trim()
            dataset_name = dataset_name ? dataset_name : 'dataset'
            dataset_name = dataset_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            def sample_name = row['id'] == null ? 'sample' : row['id'].toString().trim()
            sample_name = sample_name ? sample_name : 'sample'
            sample_name = sample_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            [
                dataset_name,
                sample_name,
                file(row['dnam_beta_matrix_file'], checkIfExists: true)
            ]
        }

    //
    // MODULE: Preprocess methylation array data with minfi
    //
    ch_minfi_samples = ch_samplesheet
        .flatMap { rows -> (rows instanceof List) ? rows : [rows] }
        .filter { row -> row['sentrix_id'] || row['sentrix_position'] || row['idats_basename'] }
        .map { row ->
            def dataset_name = row['dataset'] == null ? 'dataset' : row['dataset'].toString().trim()
            dataset_name = dataset_name ? dataset_name : 'dataset'
            dataset_name = dataset_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            def sample_name = row['id'] == null ? 'sample' : row['id'].toString().trim()
            sample_name = sample_name ? sample_name : 'sample'
            sample_name = sample_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            [
                dataset_name,
                sample_name,
                row['sentrix_id'],
                row['sentrix_position'],
                row['idats_basename']
            ]
        }

    PREPROCESS_MINFI (
        ch_minfi_samples
    )
    ch_versions = ch_versions.mix(PREPROCESS_MINFI.out.versions)

    ch_minfi_pairs = PREPROCESS_MINFI.out.corrected_inputs

    ch_pairs_for_correction = ch_minfi_pairs.mix(ch_precomputed_pairs)

    //
    // MODULE: Replace beta values with NaN where detection p-value >= threshold
    // Run per sample, then concatenate corrected betas by dataset.
    //
    P_VAL_CORRECTION (
        ch_pairs_for_correction
    )
    ch_versions = ch_versions.mix(P_VAL_CORRECTION.out.versions)

    ch_corrected_or_passthrough_betas = P_VAL_CORRECTION.out.corrected_betas
        .mix(ch_precomputed_beta_only)

    //
    // MODULE: Filter DNAm betas per sample (common probes -> missingness -> variance)
    //
    ch_common_probes = channel.value(file(params.common_probes, checkIfExists: true))

    FILTER_COMMON_PROBES (
        ch_corrected_or_passthrough_betas,
        ch_common_probes
    )
    ch_versions = ch_versions.mix(FILTER_COMMON_PROBES.out.versions)

    ch_common_filtered_by_dataset = FILTER_COMMON_PROBES.out.filtered_betas
        .groupTuple()
        .map { dataset_name, sample_names, beta_paths ->
            [
                dataset_name,
                sample_names,
                beta_paths.collect { beta_path -> beta_path.toString() }
            ]
        }

    CONCATENATE_DNAM (
        ch_common_filtered_by_dataset
    )
    ch_versions = ch_versions.mix(CONCATENATE_DNAM.out.versions)

    ch_dataset_betas_for_missing = CONCATENATE_DNAM.out.beta_matrix
        .map { beta_matrix ->
            def dataset_name = beta_matrix.baseName.replaceFirst(/\.beta_matrix$/, '')
            [dataset_name, dataset_name, beta_matrix]
        }

    FILTER_BY_MISSING (
        ch_dataset_betas_for_missing
    )
    ch_versions = ch_versions.mix(FILTER_BY_MISSING.out.versions)

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
