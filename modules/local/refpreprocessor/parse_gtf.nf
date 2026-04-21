process PARSE_GTF {
    tag "$gtf"
    label 'process_medium'
    publishDir "${params.annotation_cache_dir}", mode: 'copy', overwrite: true

    conda "bioconda::bioconductor-ballgown=2.38.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-ballgown:2.38.0--r44hdfd78af_0 '
        : 'biocontainers/bioconductor-ballgown:2.38.0--r44hdfd78af_0 '}"

    input:
    path gtf
    path saf
    val full_genome_name

    output:
    path "${full_genome_name}.annotations.full.csv", emit: annotations
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_file = "${full_genome_name}.annotations.full.csv"

    """
    #!/usr/bin/env Rscript

    library(ballgown)

    parse_gtf <- function(gtf_file_path, saf_file_path, output_file_path) {
        SAF <- read.table(saf_file_path, sep = "\t", header = TRUE)
        GeneLength <- rowsum(SAF\$End-SAF\$Start+1, SAF\$GeneID)
        annot <- gffRead(gtf_file_path)

        # Retrieve the attributes
        annot\$gene_id = getAttributeField(annot\$attributes,
            field = "gene_id")

        annot\$gene_name = getAttributeField(annot\$attributes,
            field = "gene_name")

        annot\$gene_biotype = getAttributeField(annot\$attributes,
            field = "gene_biotype")

        tokeep <- c("seqname", "feature", "gene_id", "gene_name", "gene_biotype")
        annot_filtered <- annot[tokeep]

        # Remove quotes
        annot_filtered[] <- lapply(annot_filtered, gsub, pattern='"', replacement='')

        # Keep only genes
        annot_filtered <- annot_filtered[annot_filtered\$feature == 'gene', ]

        # Convert gene lengths to dataframe
        gene_lengths <- as.data.frame(GeneLength)
        gene_lengths\$gene_id <- rownames(GeneLength)

        # Reorder columns to have gene_id as the first column
        gene_lengths <- gene_lengths[, c("gene_id", names(gene_lengths)[-ncol(gene_lengths)])]

        # Name the columns
        colnames(gene_lengths) <- c("gene_id", "gene_length")

        # Join DataFrames on 'gene_id' column
        full_annot <- merge(annot_filtered, gene_lengths, by = "gene_id", all = TRUE)

        # Rename columns
        colnames(full_annot) <- c("id", "chr", "feature", "name", "biotype", "length")

        # Drop feature column
        full_annot\$feature <- NULL

        # Strip trailing semicolons from the biotype column
        full_annot\$biotype <- sub(";\$", "", full_annot\$biotype)

        # Sort annotations by id
        full_annot <- full_annot[order(full_annot\$id), ]

        # Write processed metadata
        write.table(full_annot, output_file_path, col.names = TRUE, row.names = FALSE, sep = ',', quote = FALSE)
    }

    # Parse the GTF file
    parse_gtf("${gtf}", "${saf}", "${output_file}")

    # Create versions file
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    ballgown: "', packageVersion("ballgown"), '"')
    ), "versions.yml")
    """
}
