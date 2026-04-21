#!/usr/bin/env Rscript

# Minfi preprocessing script for methylation array data
# This script reads IDAT files based on a samplesheet and performs preprocessing

library(minfi)
library(data.table)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
opt <- list(samplesheet = NULL, outdir = ".", out_prefix = NULL)

i <- 1
while (i <= length(args)) {
    if (args[i] %in% c("-s", "--samplesheet")) {
        opt$samplesheet <- args[i + 1]
        i <- i + 2
    } else if (args[i] %in% c("-o", "--outdir")) {
        opt$outdir <- args[i + 1]
        i <- i + 2
    } else if (args[i] %in% c("--out-prefix")) {
        opt$out_prefix <- args[i + 1]
        i <- i + 2
    } else {
        i <- i + 1
    }
}

# Validate input
if (is.null(opt$samplesheet)) {
    stop("Samplesheet is required. Use --samplesheet to specify the path.")
}

if (!file.exists(opt$samplesheet)) {
    stop(paste("Samplesheet not found:", opt$samplesheet))
}

if (is.null(opt$out_prefix) || trimws(opt$out_prefix) == "") {
    stop("--out-prefix is required.")
}

cat("=== Minfi Preprocessing ===\n")
cat("Samplesheet:", opt$samplesheet, "\n")
cat("Output directory:", opt$outdir, "\n\n")

safe_sample_label <- function(targets, i) {
    preferred_cols <- c("sample", "sample_id", "sample_name", "Sample_Name", "sample_name")
    for (col in preferred_cols) {
        if (col %in% colnames(targets)) {
            value <- as.character(targets[[col]][i])
            if (!is.na(value) && value != "") {
                return(value)
            }
        }
    }
    as.character(i)
}

# Read the samplesheet (supports both CSV and TSV)
cat("Reading samplesheet...\n")
first_line <- readLines(opt$samplesheet, n = 1)
if (grepl("\t", first_line)) {
    targets <- read.delim(opt$samplesheet, stringsAsFactors = FALSE)
} else {
    targets <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)
}

cat("Found", nrow(targets), "samples\n\n")

if (nrow(targets) != 1) {
    stop("PREPROCESS_MINFI expects a one-sample samplesheet.")
}

# Ensure idats_basename column exists
if (!"idats_basename" %in% colnames(targets)) {
    stop("Samplesheet must contain an 'idats_basename' column with paths to IDAT files.")
}

# Normalize idats_basename values to avoid hidden mismatches due to whitespace.
targets$idats_basename <- trimws(as.character(targets$idats_basename))

# Build minfi Basename.
# Preferred input is idats_basename + sentrix_id + sentrix_position, where Basename is:
#   <idats_basename>/<sentrix_id>_<sentrix_position>
has_sentrix_cols <- all(c("sentrix_id", "sentrix_position") %in% colnames(targets))

if (has_sentrix_cols) {
    sentrix_id <- trimws(as.character(targets$sentrix_id))
    sentrix_position <- trimws(as.character(targets$sentrix_position))

    if (any(is.na(sentrix_id) | sentrix_id == "") || any(is.na(sentrix_position) | sentrix_position == "")) {
        stop("Samplesheet has empty sentrix_id/sentrix_position values. These are required for IDAT mode.")
    }

    targets$Basename <- file.path(targets$idats_basename, paste0(sentrix_id, "_", sentrix_position))
} else {
    # Backward-compatible fallback for direct-basename samplesheets.
    targets$Basename <- targets$idats_basename
}

# minfi requires unique basenames.
dup_idx <- duplicated(targets$Basename) | duplicated(targets$Basename, fromLast = TRUE)
if (any(dup_idx)) {
    dup_values <- unique(targets$Basename[dup_idx])
    dup_preview <- paste(head(dup_values, 5), collapse = "\n  - ")
    stop(
        paste0(
            "Duplicate IDAT basenames detected (", length(dup_values), " unique duplicates). ",
            "Each sample must map to a unique Sentrix_ID + Sentrix_Position combination.\n",
            "Examples:\n  - ", dup_preview
        )
    )
}

# Read IDAT files into RGChannelSet
cat("Reading IDAT files...\n")

# Validate file existence first
cat("Checking for IDAT file accessibility...\n")
missing_files <- NULL
for (i in seq_len(nrow(targets))) {
    basename <- targets$Basename[i]
    grn_file <- paste0(basename, "_Grn.idat")
    red_file <- paste0(basename, "_Red.idat")
    if (!file.exists(grn_file)) {
        missing_files <- c(missing_files, grn_file)
    }
    if (!file.exists(red_file)) {
        missing_files <- c(missing_files, red_file)
    }
}

if (!is.null(missing_files)) {
    stop(paste("Missing IDAT files for", length(missing_files) / 2, "samples. Examples:\n  ",
               paste(head(missing_files, 4), collapse = "\n  ")))
}
cat("All expected IDAT files found.\n\n")

cat("Validating IDAT file parseability...\n")
bad_files <- list()
for (i in seq_len(nrow(targets))) {
    basename <- targets$Basename[i]
    sample_label <- safe_sample_label(targets, i)
    grn_file <- paste0(basename, "_Grn.idat")
    red_file <- paste0(basename, "_Red.idat")

    for (idat_file in c(grn_file, red_file)) {
        parse_error <- tryCatch({
            suppressWarnings(illuminaio::readIDAT(idat_file))
            NULL
        }, error = function(e) {
            conditionMessage(e)
        })

        if (!is.null(parse_error)) {
            bad_files[[length(bad_files) + 1]] <- list(
                row = i,
                sample = sample_label,
                basename = basename,
                file = idat_file,
                error = parse_error
            )
        }
    }
}

if (length(bad_files) > 0) {
    preview <- vapply(
        bad_files[seq_len(min(length(bad_files), 5))],
        function(x) {
            paste0(
                "row ", x$row,
                " (sample=", x$sample,
                ", basename=", x$basename,
                ")\n    file: ", x$file,
                "\n    error: ", x$error
            )
        },
        character(1)
    )
    stop(
        paste0(
            "Found ", length(bad_files), " unreadable IDAT file(s).\n",
            "Examples:\n  ", paste(preview, collapse = "\n  "),
            "\nPlease remove or replace the corrupted IDAT files and rerun."
        )
    )
}
cat("All IDAT files are parseable.\n\n")


sample_targets <- targets[1, , drop = FALSE]
sample_label <- safe_sample_label(targets, 1)
dataset_label <- if ("dataset" %in% colnames(targets)) as.character(targets$dataset[1]) else "dataset"

cat("Processing sample one-by-one via task boundary...\n\n")
cat("Sample:", sample_label, "| Dataset:", dataset_label, "\n")

rg_set <- tryCatch({
    read.metharray.exp(targets = sample_targets, verbose = TRUE)
}, error = function(e) {
    err_msg <- conditionMessage(e)
    cat("Initial read failed, evaluating retry strategy...\n")
    cat("Error message:", err_msg, "\n")

    if (grepl("different array size", err_msg, ignore.case = TRUE)) {
        cat("Detected mixed array size input. Retrying with force=TRUE...\n")
        return(read.metharray.exp(targets = sample_targets, verbose = TRUE, force = TRUE))
    }

    stop(
        paste0(
            "Failed to read IDAT files for sample '", sample_label, "'. ",
            "Original error: ", err_msg
        )
    )
})

cat("Created RGChannelSet with", ncol(rg_set), "sample\n")

detP <- detectionP(rg_set)
colnames(detP) <- pData(rg_set)$sample
detp_out <- file.path(opt$outdir, paste0(opt$out_prefix, ".detection_pvalues.tsv"))
fwrite(as.data.table(detP, keep.rownames = "probe_id"), file = detp_out, sep = "\t")
cat("Saved detection p-values to", basename(detp_out), "\n")

rm(detP)
invisible(gc())

if (ncol(rg_set) < 2) {
    cat("Single-sample input detected. Using preprocessNoob fallback instead of preprocessFunnorm.\n")
    m_set <- preprocessNoob(rg_set)
} else {
    m_set <- tryCatch({
        preprocessFunnorm(rg_set)
    }, error = function(e) {
        err_msg <- conditionMessage(e)
        if (grepl("initial centers are not distinct", err_msg, ignore.case = TRUE)) {
            cat("preprocessFunnorm sex-clustering failed; retrying with preprocessNoob fallback.\n")
            return(preprocessNoob(rg_set))
        }
        stop(e)
    })
}
beta <- getBeta(m_set)
colnames(beta) <- pData(m_set)$sample
beta_out <- file.path(opt$outdir, paste0(opt$out_prefix, ".normalized_betas.tsv"))
fwrite(as.data.table(beta, keep.rownames = "probe_id"), file = beta_out, sep = "\t")
cat("Saved normalized betas to", basename(beta_out), "\n\n")

rm(m_set, rg_set, beta)
invisible(gc())

cat("=== Preprocessing complete ===\n")
