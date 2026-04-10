process CLASS_REPORTER {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/numpy_pandas_scikit-learn_seaborn_pruned:8e3b6a660fe423fc' :
        'community.wave.seqera.io/library/numpy_pandas_scikit-learn_seaborn_pruned:e21127d5908bdc2b' }"

    input:
    path subtypes_file

    output:
    path "class_distribution.svg", emit: class_distribution_svg
    path "class_distribution.csv", emit: class_distribution_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import os
    import pandas as pd

    # Fix numba caching
    os.environ['NUMBA_CACHE_DIR'] = os.path.join(os.getcwd(), '.numba')
    os.makedirs(os.environ['NUMBA_CACHE_DIR'], exist_ok=True)

    # Fix matplotlib cache directory
    os.environ['MPLCONFIGDIR'] = os.path.join(os.getcwd(), '.matplotlib')
    os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)

    import matplotlib.pyplot as plt

    def plot_class_distribution(y_series, title):
         # Use seaborn color palette
        colormap = plt.colormaps.get_cmap('tab20')
        counts = y_series.value_counts().sort_values(ascending=False)
        plt.figure(figsize=(10, 5))
        counts.plot(kind='bar', color=colormap.colors)

        for i, count in enumerate(counts):
            plt.text(i, count + 0.5, str(count), ha='center', va='bottom', rotation=90)

        plt.title(title)
        plt.xlabel("Class")
        plt.ylabel("Count")
        plt.xticks(rotation=90)
        plt.savefig("class_distribution.svg", bbox_inches='tight')

    # Read subtypes file
    subtypes_df = pd.read_csv("${subtypes_file}", sep="\t")
    plot_class_distribution(subtypes_df.iloc[:, 1], "Class Distribution")

    # Save the distributions as CSV
    subtypes_df.iloc[:, 1].value_counts().sort_values(ascending=False).to_csv("class_distribution.csv", header=["Count"])

    # Create versions file
    import pandas
    import numpy
    import seaborn
    import matplotlib
    import sklearn
    from umap import __version__ as umap_version
    import sys

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
        f.write(f'    numpy: "{numpy.__version__}"\\n')
        f.write(f'    seaborn: "{seaborn.__version__}"\\n')
        f.write(f'    matplotlib: "{matplotlib.__version__}"\\n')
        f.write(f'    sklearn: "{sklearn.__version__}"\\n')
        f.write(f'    umap: "{umap_version}"\\n')
    """
}
