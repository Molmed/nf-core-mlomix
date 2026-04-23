process FINALIZE {
    label 'process_single'

    conda "conda-forge::pandas=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.3' :
        'biocontainers/pandas:2.2.3' }"

    input:
    path classes_file
    path dnam_transposed
    path gex_transposed

    output:
    path "labels.csv", emit: labels_csv
    path "features.dnam.csv", emit: dnam_features_csv
    path "features.gex.csv", emit: gex_features_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    # Read labels (classes)
    labels = pd.read_csv("${classes_file}", sep="\t", index_col=0)
    labels = labels.iloc[:, 0]
    label_samples = labels.index.tolist()

    # Write labels as labels.csv
    labels.to_csv("labels.csv", sep=",")

    # Read and reindex DNAM features
    dnam = pd.read_csv("${dnam_transposed}", sep=",", index_col=0)
    dnam = dnam.loc[label_samples]
    dnam.to_csv("features.dnam.csv", sep=",")

    # Read and reindex GEX features (only if file provided and not empty)
    try:
        gex = pd.read_csv("${gex_transposed}", sep=",", index_col=0)
        gex = gex.loc[label_samples]
        gex.to_csv("features.gex.csv", sep=",")
    except FileNotFoundError:
        # Create empty GEX file if not provided
        pd.DataFrame().to_csv("features.gex.csv", sep=",")

    # Create versions file
    import pandas
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
    """

    stub:
    """
    touch labels.csv
    touch features.dnam.csv
    touch features.gex.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        pandas: "stub"
    END_VERSIONS
    """
}
