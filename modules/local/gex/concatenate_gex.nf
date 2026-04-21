process CONCATENATE_GEX {
    label 'process_medium'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(dataset_name), val(sample_names), path(sample_paths)

    output:
    path "${dataset_name}.csv", emit: concatenated_gex_csv
    val dataset_name, emit: dataset_name
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def file_to_sample_map = [sample_paths, sample_names].transpose().collectEntries { path, name ->
        [path.toString(), name.toString()]
    }
    def mapping_entries = file_to_sample_map.collect { k, v -> "\"${k}\": \"${v}\"" }.join(", ")
    """
    #!/usr/bin/env python3

    import pandas as pd
    import glob
    import sys
    import os
    import json

    def concatenate():
        # File to sample name mapping
        file_to_sample = {${mapping_entries}}

        # Get all input files in the work directory
        input_files = glob.glob("*.*")

        dfs = {}

        for file_path in input_files:
            file_basename = os.path.basename(file_path)

            # Get the sample name from the mapping
            if file_basename in file_to_sample:
                sample_name = file_to_sample[file_basename]
            else:
                print(f"Warning: Could not find sample name for file {file_path}")
                continue

            # Read the file
            sample_df = pd.read_csv(file_path,
                                    delim_whitespace=True,
                                    header=None,
                                    index_col=0)

            # Remove index name
            sample_df.index.name = None

            # Call the other column "count"
            sample_df.columns = ['count']

            # Create a new column with the sample name
            dfs[sample_name] = sample_df['count']

        # Concatenate all dataframes, using key as column name
        data = pd.concat(dfs, axis=1)

        # Ensure a strict numeric matrix and replace missing genes-per-sample with 0 counts.
        data = data.apply(pd.to_numeric, errors='coerce').fillna(0)

        # Sort columns by name
        data = data.sort_index(axis=1)

        # Write to file
        data.to_csv("${dataset_name}.csv")

    # Filter the annotations
    concatenate()

    # Create versions file
    import pandas
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
    """
}
