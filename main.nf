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
        random_seed
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
    PIPELINE_INITIALISATION.out.run_gex
        .subscribe { params.run_gex = it }

    PIPELINE_INITIALISATION.out.run_dnam
        .subscribe { params.run_dnam = it }

    PIPELINE_INITIALISATION.out.use_precomputed_dnam
        .subscribe { params.use_precomputed_dnam = it }

    PIPELINE_INITIALISATION.out.precomputed_beta
        .subscribe { params.precomputed_dnam_beta_matrix = it }

    PIPELINE_INITIALISATION.out.precomputed_pvals
        .subscribe { params.precomputed_dnam_pvals = it }

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
        PIPELINE_INITIALISATION.out.random_seed
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
