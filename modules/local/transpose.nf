process TRANSPOSE {
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path matrix_csv
    val suffix

    output:
    path "transposed_${suffix}.csv", emit: transposed_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    def transpose(file_path, output_suffix):
        # Load the data
        data = pd.read_csv(file_path, sep=",", index_col=0)

        # Format
        data = data.T
        data.index.name = 'id'
        data.columns.name = None

        # Dump to file
        data.to_csv(f"transposed_{output_suffix}.csv")

    # Transpose the matrix
    transpose("${matrix_csv}", "${suffix}")

    # Create versions file
    import pandas
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
    """
}
