process MERGE_DATASETS {
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path(filtered_data_files)

    output:
    path "merged.csv", emit: merged_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    def merge(files_comma_separated):
        # Split the input into a list of dataset names
        datasets = files_comma_separated.split(',')
        gex_dataframes = []
        pheno_dataframes = []

        for dataset in datasets:
            df = pd.read_csv(dataset, index_col=0)
            gex_dataframes.append(df)

        # Concatenate all DataFrames
        merged_gex_df = pd.concat(gex_dataframes, axis=1, join='inner')
        merged_gex_df.to_csv("merged.csv", index=True)

    # Filter the annotations
    merge("${filtered_data_files.join(',')}")

    # Create versions file
    import pandas
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
    """
}
