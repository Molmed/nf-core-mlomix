process CONCATENATE_DNAM {
    label 'process_medium'

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

beta_inputs = [${beta_mapping_entries}]

def extract_sample_series(file_path, sample_name):
    dataframe = pd.read_csv(file_path, sep="\t")

    if "probe_id" in dataframe.columns:
        dataframe = dataframe.set_index("probe_id")
    else:
        dataframe = dataframe.set_index(dataframe.columns[0])

    if sample_name in dataframe.columns:
        series = dataframe[sample_name]
    else:
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
    data = data.sort_index(axis=1)
    data.to_csv(output_file, sep="\t", index_label="probe_id")
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
