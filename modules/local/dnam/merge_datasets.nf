process MERGE_DATASETS {
    label 'process_high'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path beta_matrices

    output:
    path "merged.beta_matrix.tsv", emit: beta_matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

import pandas as pd
import sys
import gc
inputs = [
    path for path in [${beta_matrices.collect { matrix_path -> "\"${matrix_path.toString()}\"" }.join(', ')}]
    if path
]

if not inputs:
    raise ValueError("No dataset beta matrices were provided to MERGE_DATASETS")

merged = None
value_col_counter = 0

for input_path in inputs:
    frame = pd.read_csv(input_path, sep="\t", header=0, index_col=0)
    if frame.shape[1] < 1:
        raise ValueError(
            f"Expected at least 2 columns in {input_path}, got {frame.shape[1] + 1}"
        )

    value_col_counter += frame.shape[1]

    if merged is None:
        merged = frame
    else:
        merged = merged.join(frame, how="outer")

    del frame
    gc.collect()

if merged is None:
    raise ValueError("No data available after reading dataset beta matrices")

merged.index.name = None
merged.to_csv("merged.beta_matrix.tsv", sep="\t", header=True)

import pandas
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'    python: "{sys.version.split()[0]}"\\n')
    f.write(f'    pandas: "{pandas.__version__}"\\n')
    """

    stub:
    """
    touch merged.beta_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        pandas: "stub"
    END_VERSIONS
    """
}
