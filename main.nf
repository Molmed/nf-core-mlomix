#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/mlomix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/mlomix
    Website: https://nf-co.re/mlomix
    Slack  : https://nfcore.slack.com/channels/mlomix
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MLOMIX  } from './workflows/mlomix'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_mlomix_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_mlomix_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_mlomix_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')

def inferRunFlags(String inputFile) {
    def lines = new File(inputFile).readLines().findAll { line -> line?.trim() }
    if (!lines) {
        return [run_gex: false, run_dnam: false]
    }

    def headers = lines[0].split(',', -1).collect { header -> header.trim() }
    def run_gex = false
    def run_dnam = false

    lines.drop(1).each { line ->
        def values = line.split(',', -1)
        def row = [:]

        headers.eachWithIndex { header, index ->
            row[header] = index < values.size() ? values[index].trim() : ''
        }

        if (row.gex_feature_counts_file) {
            run_gex = true
        }
        if (row.dnam_beta_matrix_file || row.dnam_pvals_file || row.sentrix_id || row.sentrix_position || row.idats_dir) {
            run_dnam = true
        }
    }

    return [run_gex: run_gex, run_dnam: run_dnam]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_MLOMIX {

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

    //
    // WORKFLOW: Run pipeline
    //
    MLOMIX (
        ch_gex_samplesheet,
        ch_datasets,
        ch_batches,
        ch_classes,
        ch_dnam_samplesheet,
        ch_dnam_beta_matrix,
        ch_dnam_pvals,
        genome,
        annotation_version,
        random_seed,
        run_gex,
        run_dnam
    )
    emit:
    versions = MLOMIX.out.versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    def mode_info = [
        run_gex: false,
        run_dnam: false,
    ]

    if (params.input) {
        mode_info = inferRunFlags(params.input)
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_MLOMIX (
        PIPELINE_INITIALISATION.out.gex_samplesheet,
        PIPELINE_INITIALISATION.out.datasets,
        PIPELINE_INITIALISATION.out.batches,
        PIPELINE_INITIALISATION.out.classes,
        PIPELINE_INITIALISATION.out.dnam_samplesheet,
        PIPELINE_INITIALISATION.out.dnam_beta_matrix,
        PIPELINE_INITIALISATION.out.dnam_pvals,
        params.genome,
        PIPELINE_INITIALISATION.out.annotation_version,
        PIPELINE_INITIALISATION.out.random_seed,
        mode_info.run_gex,
        mode_info.run_dnam
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        channel.empty()
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
