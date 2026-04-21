process UMAP {
    label 'process_single'
    scratch true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_glasbey_matplotlib_pandas_pruned:b7aa57ab40ac88be' :
        'community.wave.seqera.io/library/pip_glasbey_matplotlib_pandas_pruned:7899227bf8a60e88' }"

    input:
    val prefix
    path gex_path
    path labels_file
    val transformed
    val random_seed

    output:
    path "${prefix}.umap.svg", emit: umap_svg
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import os
    import tempfile

    # Fix numba caching
    os.environ['NUMBA_CACHE_DIR'] = os.path.join(os.getcwd(), '.numba')
    os.makedirs(os.environ['NUMBA_CACHE_DIR'], exist_ok=True)

    # Fix matplotlib cache directory
    os.environ['MPLCONFIGDIR'] = os.path.join(os.getcwd(), '.matplotlib')
    os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)

    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler
    import umap.umap_ as umap
    import glasbey

    class GexUmap():
        def __init__(self,
                    filename,
                    counts_file,
                    transformed=False):
            self._filename = filename
            self._counts_file = counts_file
            self._FONT_SIZE = 15
            self._FIG_SIZE = (8, 8)

            # Read labels from tsv
            self._labels = pd.read_csv("${labels_file}", sep="\t", index_col=0)
            self._labels = self._labels.iloc[:, 0]
            self._transformed = transformed

        def run(self):
            # Load your dataset
            data = pd.read_csv(self._counts_file, index_col=0)

            if not self._transformed:
                data = data.T

            labels = self._labels

            # Align data and labels by common sample names
            common_samples = data.index.intersection(labels.index)
            if len(common_samples) == 0:
                raise ValueError(f"No common samples between data ({len(data)} rows) and labels ({len(labels)} rows)")

            data = data.loc[common_samples].sort_index()
            labels = labels.loc[common_samples].sort_index()

            palette = glasbey.create_palette(palette_size=len(labels.unique()))
            colormap = dict(zip(labels.unique(), palette))

            # Normalize and scale the data
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(data)

            # Apply UMAP for dimensionality reduction
            reducer = umap.UMAP(n_neighbors=15,
                                min_dist=0.1,
                                n_components=2,
                                random_state=${random_seed})
            umap_embedding = reducer.fit_transform(scaled_data)

            # Create a DataFrame for visualization
            umap_df = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'])
            umap_df['Label'] = labels.values

            # Sort labels alphabetically
            sorted_labels = sorted(labels.unique())

            # Convert labels to a categorical type with the desired order
            umap_df['Label'] = pd.Categorical(umap_df['Label'], categories=sorted_labels, ordered=True)

            # Plot UMAP
            plt.figure(figsize=(10, 8))
            sns.scatterplot(
                x='UMAP1',
                y='UMAP2',
                hue='Label',
                palette=colormap,
                data=umap_df,
                s=50
            )
            plt.title('UMAP')
            plt.legend(title='Label', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.xlabel('UMAP1')
            plt.ylabel('UMAP2')
            plt.tight_layout()

            plt.savefig(self._filename)

    transformed = True if "${transformed}" == "true" else False
    bu = GexUmap(filename="${prefix}.umap.svg",
                   counts_file="${gex_path}",
                   transformed=transformed)
    bu.run()

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
