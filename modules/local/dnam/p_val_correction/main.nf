process P_VAL_CORRECTION {
    tag "p_val_correction"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/pandas:1.5.3_cv1"

    input:
    path betas
    path detection_pvals

    output:
    path "corrected_betas.tsv"  , emit: corrected_betas
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def threshold = task.ext.threshold ?: '0.01'
    """
    p_val_correction.py \\
        --betas ${betas} \\
        --pvals ${detection_pvals} \\
        --threshold ${threshold} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch corrected_betas.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
