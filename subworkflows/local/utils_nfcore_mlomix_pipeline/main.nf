//
// Subworkflow with functionality specific to the nf-core/mlomix pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { CLASS_REPORTER            } from '../../../modules/local/mlomix/class_reporter/main'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    _input            //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/mlomix ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/mlomix/blob/main/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Shared class report, generated once before any workflow branching
    //
    ch_samplesheet_report = channel.fromPath(params.input, checkIfExists: true)
    CLASS_REPORTER(ch_samplesheet_report)
    ch_versions = ch_versions.mix(CLASS_REPORTER.out.versions)

    //
    // Create channels from combined GEX + DNAM samplesheet
    //
    def parsed_rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json").collect { row -> row[0] }
    def mode_info = validateInputSamplesheetModes(parsed_rows)

    if (mode_info.run_gex && !params.genome) {
        error("GEX input was detected but '--genome' is missing. Please provide a supported genome key (e.g. --genome GRCh38).")
    }

    if (mode_info.run_gex && !params.annotation_version) {
        error("GEX input was detected but '--annotation_version' is missing. Please provide an Ensembl release (e.g. --annotation_version 109).")
    }

    ch_rows = channel.value(parsed_rows)

    ch_rows
        .flatMap { rows -> rows }
        .filter { row -> row.gex_feature_counts_file }
        .set { ch_gex_rows }

    // Group GEX rows by dataset to build CONCATENATE_GEX input tuples.
    ch_gex_rows
        .map { row -> [ row.dataset, row ] }
        .groupTuple()
        .set { ch_gex_samplesheet }

    ch_gex_samplesheet
        .map { dataset ->
            def sample_names = dataset[1].collect { sample -> sample.id }
            def sample_paths = dataset[1].collect { sample -> file(sample.gex_feature_counts_file, checkIfExists: true) }
            [dataset[0], sample_names, sample_paths]
        }
        .set { ch_datasets }

    ch_gex_samplesheet
        .map { dataset ->
            dataset[1].collect { sample ->
                def batch_name = sample.batch ?: ''
                def suffix = batch_name ? "_${batch_name}" : ''
                "${sample.id}\t${dataset[0]}${suffix}"
            }.join('\n')
        }
        .map { content -> "sample\tbatch\n${content}\n" }
        .collectFile(name: "batches.tsv",
                     newLine: false,
                     keepHeader: true,
                     storeDir: "${params.outdir}/batch")
        .set { ch_batches }

    ch_gex_samplesheet
        .map { dataset -> dataset[1] }
        .collect()
        .map { grouped_rows ->
            def rows = grouped_rows.flatten()
            def header = "sample\tclass"
            def body = rows.collect { data -> "${data.id}\t${data['class'] ?: ''}" }.join('\n')
            "${header}\n${body}\n"
        }
        .collectFile(name: 'classes.tsv',
                     newLine: false,
                     storeDir: "${params.outdir}/class")
        .set { ch_classes }

    ch_rows
        .map { _ignored -> file(params.input, checkIfExists: true) }
        .set { ch_dnam_samplesheet }

    ch_rows
        .map { _ignored -> mode_info.precomputed_beta ? file(mode_info.precomputed_beta, checkIfExists: true) : null }
        .filter { item -> item != null }
        .set { ch_dnam_beta_matrix }

    ch_rows
        .map { _ignored -> mode_info.precomputed_pvals ? file(mode_info.precomputed_pvals, checkIfExists: true) : null }
        .filter { item -> item != null }
        .set { ch_dnam_pvals }

    ch_annotation_version = channel.value(params.annotation_version)
    ch_random_seed = channel.value(params.random_seed)

    emit:
    gex_samplesheet       = ch_gex_samplesheet
    datasets              = ch_datasets
    batches               = ch_batches
    classes               = ch_classes
    dnam_samplesheet      = ch_dnam_samplesheet
    dnam_beta_matrix      = ch_dnam_beta_matrix
    dnam_pvals            = ch_dnam_pvals
    annotation_version    = ch_annotation_version
    random_seed           = ch_random_seed
    versions              = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate modes and optional columns in combined input samplesheet
//
def validateInputSamplesheetModes(rows) {
    if (!rows || rows.isEmpty()) {
        error('Input samplesheet is empty after parsing.')
    }

    def run_gex = rows.any { row -> row.gex_feature_counts_file }

    // A sample can carry both modalities; DNAM detection must be based on DNAM columns, not on absence of GEX.
    def dnam_rows = rows.findAll { row ->
        row.dnam_beta_matrix_file || row.dnam_pvals_file || row.sentrix_id || row.sentrix_position || row.idats_dir
    }
    def run_dnam = !dnam_rows.isEmpty()

    // Validate that each DNAM sample uses exactly one mode (per-sample validation)
    dnam_rows.each { row ->
        def has_precomputed = row.dnam_beta_matrix_file || row.dnam_pvals_file
        def has_idat = row.sentrix_id || row.sentrix_position || row.idats_dir

        // Check for incomplete precomputed mode
        if ((row.dnam_beta_matrix_file && !row.dnam_pvals_file) || (!row.dnam_beta_matrix_file && row.dnam_pvals_file)) {
            error("DNAM sample '${row.id}' has incomplete precomputed mode: both dnam_beta_matrix_file and dnam_pvals_file are required together.")
        }

        // Check for incomplete IDAT mode
        if ((row.sentrix_id || row.sentrix_position || row.idats_dir) &&
            (!row.sentrix_id || !row.sentrix_position || !row.idats_dir)) {
            error("DNAM sample '${row.id}' has incomplete IDAT mode: sentrix_id, sentrix_position and idats_dir are all required together.")
        }

        // Check for mixing modes within a single sample
        if (has_precomputed && has_idat) {
            error("DNAM sample '${row.id}' mixes precomputed and IDAT modes. Each sample must use exactly one mode.")
        }

        // Check that sample has at least one DNAM mode
        if (!has_precomputed && !has_idat) {
            error("DNAM sample '${row.id}' has no valid DNAM input. Provide either (dnam_beta_matrix_file + dnam_pvals_file) or (sentrix_id + sentrix_position + idats_dir).")
        }
    }

    // Mixed-mode support: allow different samples to use different DNAM modes in the same run
    def precomputed_rows = dnam_rows.findAll { row -> row.dnam_beta_matrix_file || row.dnam_pvals_file }

    def use_precomputed_dnam = !precomputed_rows.isEmpty()
    def precomputed_beta = null
    def precomputed_pvals = null

    if (!precomputed_rows.isEmpty()) {
        // Expect single shared precomputed matrices across all precomputed DNAM rows
        def unique_beta = precomputed_rows.collect { row -> row.dnam_beta_matrix_file }.unique()
        def unique_pvals = precomputed_rows.collect { row -> row.dnam_pvals_file }.unique()
        if (unique_beta.size() != 1 || unique_pvals.size() != 1) {
            error('All DNAM precomputed rows must reference the same dnam_beta_matrix_file and dnam_pvals_file.')
        }

        precomputed_beta = unique_beta[0]
        precomputed_pvals = unique_pvals[0]
    }

    if (!run_gex && !run_dnam) {
        error('No runnable rows found. Provide gex_feature_counts_file and/or valid DNAM input columns.')
    }

    [
        run_gex: run_gex,
        run_dnam: run_dnam,
        use_precomputed_dnam: use_precomputed_dnam,
        precomputed_beta: precomputed_beta,
        precomputed_pvals: precomputed_pvals
    ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
