#!/usr/bin/env Rscript

# Minfi preprocessing script for methylation array data
# This script reads IDAT files based on a samplesheet and performs preprocessing

library(minfi)
library(data.table)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
opt <- list(samplesheet = NULL, outdir = ".")

i <- 1
while (i <= length(args)) {
    if (args[i] %in% c("-s", "--samplesheet")) {
        opt$samplesheet <- args[i + 1]
        i <- i + 2
    } else if (args[i] %in% c("-o", "--outdir")) {
        opt$outdir <- args[i + 1]
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

cat("=== Minfi Preprocessing ===\n")
cat("Samplesheet:", opt$samplesheet, "\n")
cat("Output directory:", opt$outdir, "\n\n")

# Read the samplesheet (supports both CSV and TSV)
cat("Reading samplesheet...\n")
first_line <- readLines(opt$samplesheet, n = 1)
if (grepl("\t", first_line)) {
    targets <- read.delim(opt$samplesheet, stringsAsFactors = FALSE)
} else {
    targets <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)
}

cat("Found", nrow(targets), "samples\n\n")

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

# Use SerialParam for sequential reading to avoid BiocParallel race conditions
# This is slower but more robust for large datasets with potential file read issues
library(BiocParallel)

result <- tryCatch({
    # Attempt parallel read first (faster for small-to-medium datasets)
    cat("Attempting parallel read...\n")
    read.metharray.exp(targets = targets, verbose = TRUE)
}, error = function(e) {
    # Fall back to serial reading on error
    cat("Parallel read failed, switching to serial mode...\n")
    cat("Error message:", conditionMessage(e), "\n")
    read.metharray.exp(targets = targets, verbose = TRUE,
                       BPPARAM = BiocParallel::SerialParam())
})

rg_set <- result
cat("Created RGChannelSet with", ncol(rg_set), "samples\n\n")

# Save RGChannelSet
saveRDS(rg_set, file = file.path(opt$outdir, "rgset.rds"))
cat("Saved RGChannelSet to rgset.rds\n")

# Generate QC report
cat("Generating QC report...\n")
pdf(file.path(opt$outdir, "qc_report.pdf"), width = 10, height = 8)

# QC plot
qc <- getQC(preprocessRaw(rg_set))
plotQC(qc)

# Density plot of beta values (raw)
densityPlot(rg_set, main = "Raw Beta Values Density")

dev.off()
cat("Saved QC report to qc_report.pdf\n\n")

# Calculate detection p-values
cat("Computing detection p-values...\n")
detP <- detectionP(rg_set)
colnames(detP) <- pData(rg_set)$sample
cat("Detection p-values matrix:", nrow(detP), "probes x", ncol(detP), "samples\n")

# Save the detection p-values to separate file
fwrite(as.data.table(detP, keep.rownames = "probe_id"), file = file.path(opt$outdir, "detection_pvalues.tsv"), sep = "\t")
cat("Saved detection p-values to detection_pvalues.tsv\n\n")

# Preprocess with functional normalization (includes Noob + dye bias correction)
# This returns GenomicRatioSet and therefore ratioConvert not needed
cat("Performing functional normalization...\n")
m_set <- preprocessFunnorm(rg_set)

# Extract beta values that are baked in the object (GenomicRatioSet) produced in preprocessFunnorm()
beta <- getBeta(m_set)
colnames(beta) <- pData(m_set)$sample
cat("Extracted beta values matrix:", nrow(beta), "probes x", ncol(beta), "samples\n")

# Save normalized beta values
fwrite(as.data.table(beta, keep.rownames = "probe_id"), file = file.path(opt$outdir, "normalized_betas.tsv"), sep = "\t")
cat("Saved normalized beta values to normalized_betas.tsv\n")

cat("\n=== Preprocessing complete ===\n")
