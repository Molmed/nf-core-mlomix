process TSNE {
    label 'process_single'
    scratch true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_glasbey_matplotlib_pandas_pruned:b7aa57ab40ac88be' :
        'community.wave.seqera.io/library/pip_glasbey_matplotlib_pandas_pruned:7899227bf8a60e88' }"

    input:
    val prefix
    path features_path
    path labels_file
    val transformed
    val random_seed

    output:
    path "${prefix}.tsne.svg", emit: tsne_svg
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
    from sklearn.manifold import TSNE as SkTSNE
    import glasbey

    class FeaturesTsne():
        def __init__(self,
                    filename,
                    features_file,
                    transformed=False):
            self._filename = filename
            self._features_file = features_file
            self._FONT_SIZE = 15
            self._FIG_SIZE = (8, 8)

            # Read labels from tsv
            self._labels = pd.read_csv("${labels_file}", sep="\t", index_col=0)
            self._labels = self._labels.iloc[:, 0]
            self._transformed = transformed

        def run(self):
            # Load feature matrix
            data = pd.read_csv(self._features_file, index_col=0)

            if not self._transformed:
                data = data.T

            labels = self._labels

            # Align data and labels by common sample names
            common_samples = data.index.intersection(labels.index)
            if len(common_samples) == 0:
                raise ValueError(f"No common samples between data ({len(data)} rows) and labels ({len(labels)} rows)")

            data = data.loc[common_samples].sort_index()
            labels = labels.loc[common_samples].sort_index()

            # DNAm matrices can contain missing values after filtering; coerce and impute for t-SNE.
            data = data.apply(pd.to_numeric, errors='coerce')
            data = data.replace([np.inf, -np.inf], np.nan)

            # Drop samples/probes that are completely missing, then impute remaining missing values.
            data = data.dropna(axis=0, how='all')
            labels = labels.loc[data.index]
            data = data.dropna(axis=1, how='all')

            if data.shape[0] < 2:
                raise ValueError("Need at least 2 samples with non-missing values for t-SNE")
            if data.shape[1] < 1:
                raise ValueError("No usable features remain after removing all-missing columns")

            data = data.fillna(data.median(axis=0))
            data = data.fillna(0)

            palette = glasbey.create_palette(palette_size=len(labels.unique()))
            colormap = dict(zip(labels.unique(), palette))

            # Normalize and scale the data
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(data)

            # Apply t-SNE for dimensionality reduction
            reducer = SkTSNE(n_components=2,
                             perplexity=30,
                             learning_rate='auto',
                             init='pca',
                             random_state=${random_seed})
            tsne_embedding = reducer.fit_transform(scaled_data)

            # Create a DataFrame for visualization
            tsne_df = pd.DataFrame(tsne_embedding, columns=['tSNE1', 'tSNE2'])
            tsne_df['Label'] = labels.values

            # Sort labels alphabetically
            sorted_labels = sorted(labels.unique())

            # Convert labels to a categorical type with the desired order
            tsne_df['Label'] = pd.Categorical(tsne_df['Label'], categories=sorted_labels, ordered=True)

            # Plot t-SNE
            plt.figure(figsize=(10, 8))
            sns.scatterplot(
                x='tSNE1',
                y='tSNE2',
                hue='Label',
                palette=colormap,
                data=tsne_df,
                s=50
            )
            plt.title('t-SNE')
            plt.legend(title='Label', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.xlabel('tSNE1')
            plt.ylabel('tSNE2')
            plt.tight_layout()

            plt.savefig(self._filename)

    transformed = True if "${transformed}" == "true" else False
    bu = FeaturesTsne(filename="${prefix}.tsne.svg",
                      features_file="${features_path}",
                      transformed=transformed)
    bu.run()

    # Create versions file
    import pandas
    import numpy
    import seaborn
    import matplotlib
    import sklearn
    import sys

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
        f.write(f'    numpy: "{numpy.__version__}"\\n')
        f.write(f'    seaborn: "{seaborn.__version__}"\\n')
        f.write(f'    matplotlib: "{matplotlib.__version__}"\\n')
        f.write(f'    sklearn: "{sklearn.__version__}"\\n')
    """
}
