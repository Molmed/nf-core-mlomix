process FILTER_BY_COMMON_MISSING_AND_VARIANCE {
    tag "filter_by_common_missing_and_variance"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/matplotlib:3.7.1"

    input:
    path betas
    path probe_list

    output:
    path "variance_filtered_betas.tsv"                , emit: variance_filtered_betas
    path "variance_distribution_before_filtering.png" , emit: variance_plot_before_png
    path "variance_distribution_before_filtering.svg" , emit: variance_plot_before_svg
    path "variance_distribution_after_filtering.png"  , emit: variance_plot_after_png
    path "variance_distribution_after_filtering.svg"  , emit: variance_plot_after_svg
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    filter_by_common_missing_and_variance.py \
        --betas ${betas} \
        --probes ${probe_list} \
        --missing_threshold ${params.missing_threshold} \
        --variance_threshold ${params.variance_threshold} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch variance_filtered_betas.tsv
    touch variance_distribution_before_filtering.png
    touch variance_distribution_before_filtering.svg
    touch variance_distribution_after_filtering.png
    touch variance_distribution_after_filtering.svg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
