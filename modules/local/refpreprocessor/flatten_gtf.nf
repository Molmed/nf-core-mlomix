process FLATTEN_GTF {
    tag "$gtf"
    label 'process_medium'

    conda "bioconda::subread=2.1.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0'
        : 'biocontainers/subread:2.1.1--h577a1d6_0'}"

    input:
    path gtf
    val full_genome_name

    output:
    path "${full_genome_name}.saf", emit: saf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output_file = "${full_genome_name}.saf"

    """
    flattenGTF ${args} -a ${gtf} -o ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(flattenGTF 2>&1 | grep "flattenGTF Version" | sed -n 's/.*Version \\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """
}
