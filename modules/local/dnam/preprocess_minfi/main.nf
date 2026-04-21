process PREPROCESS_MINFI {
    tag "minfi_preprocessing"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicanno.ilm10b4.hg19_bioconductor-illuminahumanmethylationepicmanifest_pruned:a2c54946eaf79f9e' :
        'community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicanno.ilm10b4.hg19_bioconductor-illuminahumanmethylationepicmanifest_pruned:98e7c81d1a064d77' }"

    input:
    path samplesheet

    output:
    path "*.normalized_betas.tsv"  , emit: betas
    path "*.detection_pvalues.tsv" , emit: detection_pvals
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    preprocess_minfi.R \\
        --samplesheet ${samplesheet} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
        bioconductor-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch sample_001.normalized_betas.tsv
    touch sample_001.detection_pvalues.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
        bioconductor-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """
}
