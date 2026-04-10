process DOWNLOAD_GTF {
    cache 'lenient'
    label "process_low"

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
        val genome
        val annotation_version

    output:
        path "Homo_sapiens.${genome}.${annotation_version}.gtf.gz", emit: gtf_file
        path 'versions.yml', emit: versions

    script:
    def filename="Homo_sapiens.${genome}.${annotation_version}.gtf.gz"
    def url="http://ftp.ensembl.org/pub/release-${annotation_version}/gtf/homo_sapiens/${filename}"
    """
    wget $url -O $filename
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3)
    END_VERSIONS
    """
}
