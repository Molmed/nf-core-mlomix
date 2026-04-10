process TRANSPOSE {
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path normalized_csv

    output:
    path "transposed.csv", emit: transposed_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    def transpose(file_path):
        # Load the data
        data = pd.read_csv(file_path, index_col=0)

        # Format
        data = data.T
        data.index.name = 'id'
        data.columns.name = None

        # Dump to file
        data.to_csv("transposed.csv")

    # Filter the annotations
    transpose("${normalized_csv}")

    # Create versions file
    import pandas
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
    """
}
