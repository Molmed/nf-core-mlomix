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

    batch_correct <- function(gex_path, sample_names, batch_df, subtypes, output_path) {
        x <- read.csv(gex_path, row.names = 1, header= TRUE, check.names = FALSE)

        # batch_df must contain sample and batch columns.
        # The batch value already includes dataset context from pipeline initialization.
        if (!("sample" %in% colnames(batch_df)) || !("batch" %in% colnames(batch_df))) {
            stop("batches_file must contain 'sample' and 'batch' columns")
        }

        # Ensure batch vector is ordered to match 'x'
        ordered_idx <- match(colnames(x), batch_df\$sample)
        if (any(is.na(ordered_idx))) {
            stop("Some expression matrix columns are missing in batches_file sample column")
        }
        batch <- as.character(batch_df\$batch[ordered_idx])

        # Do batch correction only when there is more than one combined batch
        # and each combined batch has at least 2 samples in the matched expression set.
        has_batch_variation <- length(unique(batch)) > 1 && all(table(batch) >= 2)

        if (!has_batch_variation) {
            write.csv(x, output_path)
            return()
        }

        # ComBat_seq expects a non-negative integer count matrix without NA/Inf.
        counts <- as.matrix(x)
        storage.mode(counts) <- "numeric"

        non_finite_mask <- !is.finite(counts)
        non_finite_n <- sum(non_finite_mask)
        if (non_finite_n > 0) {
            cat("Replacing", non_finite_n, "non-finite count values (NA/Inf) with 0 before ComBat_seq.\n")
            counts[non_finite_mask] <- 0
        }

        negative_mask <- counts < 0
        negative_n <- sum(negative_mask)
        if (negative_n > 0) {
            cat("Clamping", negative_n, "negative count values to 0 before ComBat_seq.\n")
            counts[negative_mask] <- 0
        }

        counts <- round(counts)

        # Use group parameter only if matched subtypes have >1 non-empty unique value
        group <- NULL
        if (!is.null(subtypes)) {
            matched_group <- as.character(subtypes[ordered_idx])
            valid_group <- !is.na(matched_group) & matched_group != ""
            if (sum(valid_group) > 1 && length(unique(matched_group[valid_group])) > 1) {
                group <- matched_group
            }
        }

        # Call ComBat_seq with or without group parameter
        if (!is.null(group)) {
            correcteddata <- ComBat_seq(counts = counts, batch=batch, group=group)
        } else {
            correcteddata <- ComBat_seq(counts = counts, batch=batch)
        }
        write.csv(correcteddata, output_path)
    }

    batches_df <- read.csv("${batches_file}", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    sample_names <- batches_df\$sample

    subtypes_df <- read.csv("${subtypes_file}", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    subtypes <- subtypes_df\$class

    batch_correct("${gex_path}", sample_names, batches_df, subtypes, "batch_corrected.csv")

    # Create versions file
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    sva: "', packageVersion("sva"), '"')
    ), "versions.yml")
    """
}
