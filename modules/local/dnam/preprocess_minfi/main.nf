process PREPROCESS_MINFI {
    tag "minfi_preprocessing"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicanno.ilm10b4.hg19_bioconductor-illuminahumanmethylationepicmanifest_pruned:a2c54946eaf79f9e' :
        'community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicanno.ilm10b4.hg19_bioconductor-illuminahumanmethylationepicmanifest_pruned:98e7c81d1a064d77' }"

    input:
    tuple val(dataset_name), val(sample_name), val(sentrix_id), val(sentrix_position), val(idats_basename)

    output:
    tuple val(dataset_name), val(sample_name), path("${dataset_name}__${sample_name}.normalized_betas.tsv"), path("${dataset_name}__${sample_name}.detection_pvalues.tsv"), emit: corrected_inputs
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample_sheet_name = "${dataset_name}__${sample_name}.samplesheet.csv"
    """
    cat <<'EOF' > ${sample_sheet_name}
sample,dataset,sentrix_id,sentrix_position,idats_basename
${sample_name},${dataset_name},${sentrix_id},${sentrix_position},${idats_basename}
EOF

    preprocess_minfi.R \\
        --samplesheet ${sample_sheet_name} \
        --out-prefix ${dataset_name}__${sample_name} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
        bioconductor-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch ${dataset_name}__${sample_name}.normalized_betas.tsv
    touch ${dataset_name}__${sample_name}.detection_pvalues.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
        bioconductor-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """
}
