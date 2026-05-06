process FILTER_BY_VARIANCE {
    tag "filter_by_variance"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/matplotlib:3.7.1"

    input:
    tuple val(dataset_name), val(sample_name), path(betas)

    output:
    tuple val(dataset_name), val(sample_name), path("${dataset_name}__${sample_name}.variance_filtered_betas.csv"), emit: variance_filtered_betas
    path "${dataset_name}__${sample_name}.variance_distribution_before_filtering.png"  , emit: variance_plot_before_png
    path "${dataset_name}__${sample_name}.variance_distribution_before_filtering.svg"  , emit: variance_plot_before_svg
    path "${dataset_name}__${sample_name}.variance_distribution_after_filtering.png"   , emit: variance_plot_after_png
    path "${dataset_name}__${sample_name}.variance_distribution_after_filtering.svg"   , emit: variance_plot_after_svg
    path "versions.yml"                                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cli_args = [
        "--betas ${betas}",
        "--variance_threshold ${params.dnam_variance_threshold}"
    ]
    if (params.dnam_variance_keep_top_sites != null) {
        cli_args << "--keep_top_sites ${params.dnam_variance_keep_top_sites}"
    }
    if (args) {
        cli_args << args
    }
    def cli_args_block = cli_args.join(' \\\n+        ')
    def beta_output_name = "${dataset_name}__${sample_name}.variance_filtered_betas.csv"
    def before_png_name = "${dataset_name}__${sample_name}.variance_distribution_before_filtering.png"
    def before_svg_name = "${dataset_name}__${sample_name}.variance_distribution_before_filtering.svg"
    def after_png_name = "${dataset_name}__${sample_name}.variance_distribution_after_filtering.png"
    def after_svg_name = "${dataset_name}__${sample_name}.variance_distribution_after_filtering.svg"
    """
    filter_by_variance.py \
        ${cli_args_block}

    mv variance_filtered_betas.csv ${beta_output_name}
    mv variance_distribution_before_filtering.png ${before_png_name}
    mv variance_distribution_before_filtering.svg ${before_svg_name}
    mv variance_distribution_after_filtering.png ${after_png_name}
    mv variance_distribution_after_filtering.svg ${after_svg_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${dataset_name}__${sample_name}.variance_filtered_betas.csv
    touch ${dataset_name}__${sample_name}.variance_distribution_before_filtering.png
    touch ${dataset_name}__${sample_name}.variance_distribution_before_filtering.svg
    touch ${dataset_name}__${sample_name}.variance_distribution_after_filtering.png
    touch ${dataset_name}__${sample_name}.variance_distribution_after_filtering.svg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
