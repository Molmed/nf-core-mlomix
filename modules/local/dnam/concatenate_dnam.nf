process CONCATENATE_DNAM {
    label 'process_high'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(dataset_name), val(sample_names), val(beta_matrix_paths)

    output:
    path "${dataset_name}.beta_matrix.tsv", emit: beta_matrix
    val dataset_name, emit: dataset_name
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def beta_mapping_entries = [sample_names, beta_matrix_paths].transpose().collect { sample_name, file_path ->
        "[\"${file_path.toString()}\", \"${sample_name.toString()}\"]"
    }.join(", ")
    """
    #!/usr/bin/env python3

import pandas as pd
import sys
import gc

beta_inputs = [${beta_mapping_entries}]

def extract_sample_series(file_path, sample_name):
    dataframe = pd.read_csv(file_path, sep="\t", header=None)

    if dataframe.shape[1] < 2:
        raise ValueError(
            f"Expected at least 2 columns in {file_path}, got {dataframe.shape[1]}"
        )

    if str(dataframe.iloc[0, 0]).strip().lower() == "probe_id":
        dataframe = dataframe.iloc[1:, :]

    dataframe.columns = ["probe_id"] + [f"value_{i}" for i in range(1, dataframe.shape[1])]
    dataframe = dataframe.set_index("probe_id")

    value_columns = list(dataframe.columns)
    if len(value_columns) != 1:
        raise ValueError(
            f"Could not determine a single value column for {file_path} and sample {sample_name}"
        )
    series = dataframe[value_columns[0]]

    series.index.name = None
    return series

def concatenate(inputs, output_file):
    if not inputs:
        return False

    sample_frames = {}

    for file_path, sample_name in inputs:
        sample_frames[sample_name] = extract_sample_series(file_path, sample_name)

    data = pd.concat(sample_frames, axis=1)
    data.index.name = None
    data = data.sort_index(axis=1)
    data.to_csv(output_file, sep="\t", header=False)
    return True

concatenate(beta_inputs, "${dataset_name}.beta_matrix.tsv")

import pandas
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'    python: "{sys.version.split()[0]}"\\n')
    f.write(f'    pandas: "{pandas.__version__}"\\n')
    """

    stub:
    """
    touch ${dataset_name}.beta_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        pandas: "stub"
    END_VERSIONS
    """
}

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

inputs = [
    path for path in [${beta_matrices.collect { matrix_path -> "\"${matrix_path.toString()}\"" }.join(', ')}]
    if path
]

if not inputs:
    raise ValueError("No dataset beta matrices were provided to MERGE_DATASETS")

merged = None
value_col_counter = 0

for input_path in inputs:
    frame = pd.read_csv(input_path, sep="\t", header=None, index_col=0)
    if frame.shape[1] < 1:
        raise ValueError(
            f"Expected at least 2 columns in {input_path}, got {frame.shape[1] + 1}"
        )

    frame.columns = [
        f"value_{i}"
        for i in range(value_col_counter + 1, value_col_counter + frame.shape[1] + 1)
    ]
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
merged.to_csv("merged.beta_matrix.tsv", sep="\t", header=False)

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
