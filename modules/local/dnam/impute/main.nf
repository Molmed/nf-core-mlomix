process IMPUTE {
    tag "impute"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "python:3.11-slim"

    input:
    tuple val(dataset_name), val(sample_name), path(beta_matrix)

    output:
    tuple val(dataset_name), val(sample_name), path("${dataset_name}__${sample_name}.imputed_betas.csv"), emit: imputed_betas
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import numpy as np
    from sklearn.impute import SimpleImputer
    import sklearn
    import sys

    beta_matrix = pd.read_csv("${beta_matrix}", index_col=0)
    beta_matrix = beta_matrix.apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan)

    # Median imputation (more robust to outliers)
    imp = SimpleImputer(strategy='median')
    beta_imputed = imp.fit_transform(beta_matrix)
    beta_imputed = pd.DataFrame(beta_imputed, index=beta_matrix.index, columns=beta_matrix.columns)
    beta_imputed.to_csv("${dataset_name}__${sample_name}.imputed_betas.csv")

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pd.__version__}"\\n')
        f.write(f'    numpy: "{np.__version__}"\\n')
        f.write(f'    sklearn: "{sklearn.__version__}"\\n')
    """

    stub:
    """
    touch ${dataset_name}__${sample_name}.imputed_betas.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)" )
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)" )
        sklearn: \$(python3 -c "import sklearn; print(sklearn.__version__)" )
    END_VERSIONS
    """
}
