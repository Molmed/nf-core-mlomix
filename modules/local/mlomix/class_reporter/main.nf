process CLASS_FILTER_AND_REPORT {
    tag "class_filter_and_report"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/pandas:1.5.3_cv1"

    input:
    path samplesheet

    output:
    path "class_distribution.gex.all.svg"     , emit: class_plot_gex_all
    path "class_distribution.gex.filtered.svg", emit: class_plot_gex_filtered
    path "class_distribution.dnam.all.svg"    , emit: class_plot_dnam_all
    path "class_distribution.dnam.filtered.svg", emit: class_plot_dnam_filtered
    path "classes.gex.all.tsv"                , emit: classes_gex_all
    path "classes.gex.filtered.tsv"           , emit: classes_gex_filtered
    path "classes.dnam.all.tsv"               , emit: classes_dnam_all
    path "classes.dnam.filtered.tsv"          , emit: classes_dnam_filtered
    path "samplesheet.filtered.tsv"           , emit: samplesheet_filtered
    path "kept_sample_ids.txt"                , emit: kept_sample_ids
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${projectDir}/bin/class_filter_and_report.py ${samplesheet} --gex-min ${params.gex_min_samples_per_class} --dnam-min ${params.dnam_min_samples_per_class}
    """

    stub:
    """
    touch class_distribution.gex.all.svg
    touch class_distribution.gex.filtered.svg
    touch class_distribution.dnam.all.svg
    touch class_distribution.dnam.filtered.svg
    touch classes.gex.all.tsv
    touch classes.gex.filtered.tsv
    touch classes.dnam.all.tsv
    touch classes.dnam.filtered.tsv
    touch samplesheet.filtered.tsv
    touch kept_sample_ids.txt
    echo "versions:" > versions.yml
    """
}
