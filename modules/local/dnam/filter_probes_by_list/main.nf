process FILTER_PROBES_BY_LIST {
    tag "filter_probes_by_list"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/pandas:1.5.3_cv1"

    input:
    tuple val(dataset_name), val(sample_name), path(betas)
    path probe_list

    output:
    tuple val(dataset_name), val(sample_name), path("${dataset_name}__${sample_name}.filtered_betas.tsv"), emit: filtered_betas
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output_name = "${dataset_name}__${sample_name}.filtered_betas.tsv"
    """
    filter_probes_by_list.py \
        --sample-name ${sample_name} \
        --betas ${betas} \
        --probes ${probe_list} \
        ${args}

    mv filtered_betas.tsv ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${dataset_name}__${sample_name}.filtered_betas.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
