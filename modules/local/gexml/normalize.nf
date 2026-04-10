process NORMALIZE {
    label 'process_medium'

    conda "bioconda::bioconductor-edger=4.4.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-edger:4.4.0--r44h3df3fcb_0 '
        : 'biocontainers/bioconductor-edger:4.4.0--r44h3df3fcb_0 '}"

    input:
    path merged_data_path
    path ref_path

    output:
    path "normalized.csv", emit: normalized_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(edgeR)

    # create a function `get_cpm`
    normalize <- function(gex_path, ref_path, output_path) {
        x <- read.csv(gex_path, row.names = 1, header= TRUE, check.names = FALSE)
        annot <- read.csv(ref_path, row.names = 1, header= TRUE, check.names = FALSE)

        x_length_norm <- ( (x*10^3 )/annot\$length)
        d <- DGEList(counts=x_length_norm)
        TMM <- calcNormFactors(d, method="TMM")
        CPM <- cpm(TMM, log = TRUE)
        write.csv(CPM, output_path)
    }

    normalize("${merged_data_path}", "${ref_path}", "normalized.csv")

    # Create versions file
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    E: "', packageVersion("edgeR"), '"')
    ), "versions.yml")
    """
}
