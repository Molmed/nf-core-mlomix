process BATCH_CORRECT {
    label 'process_high'

    conda "bioconda::bioconductor-sva=3.54.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-sva:3.54.0--r44h3df3fcb_1 '
        : 'biocontainers/bioconductor-sva:3.54.0--r44h3df3fcb_1 '}"

    input:
    path gex_path
    path batches_file
    path subtypes_file

    output:
    path "batch_corrected.csv", emit: batch_corrected_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(sva)

    batch_correct <- function(gex_path, sample_names, batch_names, subtypes, output_path) {
        x <- read.csv(gex_path, row.names = 1, header= TRUE, check.names = FALSE)
        # Do batch correction only if there is more than one batch and
        # more than one sample in each batch
        if (length(unique(batch_names)) <= 1 || any(table(batch_names) < 2)) {
            write.csv(x, output_path)
            return()
        }

        batches <- data.frame(batch = batch_names, row.names = sample_names)
        # Ensure 'batches' is ordered to match 'x'
        batches <- batches[colnames(x), , drop = FALSE]
        batch = batches\$batch

        # Process subtypes the same way
        subtypes <- data.frame(subtype = subtypes, row.names = sample_names)
        subtypes <- subtypes[colnames(x), , drop = FALSE]
        subtype = subtypes\$subtype

        correcteddata <- ComBat_seq(as.matrix(x), batch=batch, group=subtype)
        write.csv(correcteddata, output_path)
    }

    batches_df <- read.csv("${batches_file}", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    sample_names <- batches_df\$sample
    batches <- batches_df\$batch

    subtypes_df <- read.csv("${subtypes_file}", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    subtypes <- subtypes_df\$class

    batch_correct("${gex_path}", sample_names, batches, subtypes, "batch_corrected.csv")

    # Create versions file
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    sva: "', packageVersion("sva"), '"')
    ), "versions.yml")
    """
}
