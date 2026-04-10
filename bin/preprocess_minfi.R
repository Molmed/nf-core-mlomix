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

# Ensure idats_dir column exists
if (!"idats_dir" %in% colnames(targets)) {
    stop("Samplesheet must contain an 'idats_dir' column with paths to IDAT files.")
}

# minfi expects the IDAT path column to be named Basename.
targets$Basename <- targets$idats_dir

# Read IDAT files into RGChannelSet
cat("Reading IDAT files...\n")
rg_set <- read.metharray.exp(targets = targets)
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
