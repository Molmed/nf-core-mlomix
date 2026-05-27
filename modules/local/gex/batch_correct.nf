process BATCH_CORRECT {
    label 'process_high'

    conda "bioconda::bioconductor-sva=3.54.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-sva:3.54.0--r44h3df3fcb_1 '
        : 'biocontainers/bioconductor-sva:3.54.0--r44h3df3fcb_1 '}"

    input:
    path gex_path
    path batches_file
    path classes_file

    output:
    path "batch_corrected.csv", emit: batch_corrected_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(sva)

    cat("BATCH_CORRECT: gex_use_combat_seq_group = ${params.gex_use_combat_seq_group}\n")

    batch_correct <- function(gex_path, sample_names, batch_df, classes, use_group_str, output_path) {
        # Parse the use_group parameter properly - handle both string and boolean values from Nextflow
        use_group <- as.logical(use_group_str)
        if (is.na(use_group)) {
            use_group <- FALSE
        }

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

        # ComBat_seq cannot handle singleton batches, so keep them unchanged
        # and only correct the samples that belong to batches with >= 2 samples.
        batch_counts <- table(batch)
        singleton_batches <- names(batch_counts[batch_counts < 2])
        singleton_mask <- batch %in% singleton_batches
        corrected_mask <- !singleton_mask

        has_batch_variation <- length(unique(batch[corrected_mask])) > 1

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

        # Use group parameter only when use_group is TRUE and matched classes have >1 non-empty unique value.
        group <- NULL
        full_mod <- FALSE
        if (use_group) {
            matched_group <- as.character(classes[ordered_idx])
            valid_group <- !is.na(matched_group) & matched_group != ""
            if (sum(valid_group) > 1 && length(unique(matched_group[valid_group])) > 1) {
                group <- matched_group
                full_mod <- TRUE
            }
        }

        if (any(singleton_mask)) {
            cat("Excluding", sum(singleton_mask), "samples from ComBat_seq because their batches are singletons.\n")
        }

        counts_for_correction <- counts[, corrected_mask, drop = FALSE]
        batch_for_correction <- batch[corrected_mask]
        group_for_correction <- if (is.null(group)) NULL else group[corrected_mask]

        cat("Performing ComBat_seq batch correction with use_group =", use_group, "and full_mod =", full_mod, "\n")

        correcteddata <- ComBat_seq(
            counts = counts_for_correction,
            batch = batch_for_correction,
            group = group_for_correction,
            full_mod = full_mod
        )

        output_counts <- counts
        output_counts[, corrected_mask] <- correcteddata
        write.csv(output_counts, output_path)
    }

    batches_df <- read.csv("${batches_file}", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    sample_names <- batches_df\$sample

    classes_df <- read.csv("${classes_file}", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    classes <- classes_df\$class

    batch_correct("${gex_path}", sample_names, batches_df, classes, "${params.gex_use_combat_seq_group}", "batch_corrected.csv")

    # Create versions file
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    sva: "', packageVersion("sva"), '"')
    ), "versions.yml")
    """
}
