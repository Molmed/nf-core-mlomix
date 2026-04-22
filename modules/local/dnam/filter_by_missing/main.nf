process FILTER_BY_MISSING {
    tag "filter_by_missing"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/pandas:1.5.3_cv1"

    input:
    tuple val(dataset_name), val(sample_name), path(betas)

    output:
    tuple val(dataset_name), val(sample_name), path("${dataset_name}__${sample_name}.missing_filtered_betas.csv"), emit: missing_filtered_betas
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output_name = "${dataset_name}__${sample_name}.missing_filtered_betas.csv"
    """
    filter_by_missing.py \
        --betas ${betas} \
        --missing_threshold ${params.missing_threshold} \
        ${args}

    mv missing_filtered_betas.csv ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${dataset_name}__${sample_name}.missing_filtered_betas.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
