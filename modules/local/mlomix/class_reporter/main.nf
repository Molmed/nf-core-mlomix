process CLASS_REPORTER {
    tag "class_reporter"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/pandas:1.5.3_cv1"

    input:
    path samplesheet

    output:
    path "class_distribution.gex.svg"  , emit: class_plot_gex
    path "class_distribution.dnam.svg" , emit: class_plot_dnam
    path "classes.gex.tsv"            , emit: classes_gex
    path "classes.dnam.tsv"           , emit: classes_dnam
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import os
    import sys
    import pandas as pd

    # Fix matplotlib cache directory
    os.environ['MPLCONFIGDIR'] = os.path.join(os.getcwd(), '.matplotlib')
    os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def plot_class_distribution(y_series, title, outname):
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
        plt.savefig(outname, bbox_inches='tight')

    def plot_empty(message, title, outname):
        plt.figure(figsize=(10, 5))
        plt.text(0.5, 0.5, message, ha='center', va='center', fontsize=14, transform=plt.gca().transAxes)
        plt.title(title)
        plt.xlabel("Class")
        plt.ylabel("Count")
        plt.savefig(outname, bbox_inches='tight')

    # Auto-detect delimiter
    samplesheet_path = "${samplesheet}"
    if samplesheet_path.endswith('.tsv'):
        sep = '\t'
    else:
        sep = ','

    df = pd.read_csv(samplesheet_path, sep=sep)

    # Determine modality by presence of GEX/DNAM-specific columns
    has_gex = df['gex_feature_counts_file'].notna() if 'gex_feature_counts_file' in df.columns else pd.Series([False] * len(df))
    has_dnam = (df['dnam_beta_matrix_file'].notna() | df['dnam_pvals_file'].notna() | df['idats_basename'].notna()) if any(c in df.columns for c in ['dnam_beta_matrix_file','dnam_pvals_file','idats_basename']) else pd.Series([False] * len(df))

    # For each modality produce TSV and SVG
    def handle_modality(mask, modality_label, tsv_name, svg_name):
        subdf = df[mask]
        if 'class' not in subdf.columns:
            print(f"WARNING: Samplesheet does not contain a 'class' column for {modality_label}. Generating empty report.")
            plot_empty(f"No 'class' column found for {modality_label}", f"Class Distribution ({modality_label})", svg_name)
            pd.DataFrame(columns=["class", "Count"]).to_csv(tsv_name, sep='\t', index=False)
            return
        classes = subdf['class'].dropna()
        if classes.empty:
            print(f"WARNING: 'class' column present for {modality_label} but contains no values. Generating empty report.")
            plot_empty(f"'class' column has no values for {modality_label}", f"Class Distribution ({modality_label})", svg_name)
            pd.DataFrame(columns=["class", "Count"]).to_csv(tsv_name, sep='\t', index=False)
        else:
            plot_class_distribution(classes, f"Class Distribution ({modality_label})", svg_name)
            vc = classes.value_counts().sort_values(ascending=False)
            vc.to_csv(tsv_name, sep='\t', header=["Count"], index_label='class')

    # Handle GEX
    if has_gex.any():
        handle_modality(has_gex, 'GEX', 'classes.gex.tsv', 'class_distribution.gex.svg')
    else:
        # no GEX samples
        plot_empty('No GEX samples found', 'Class Distribution (GEX)', 'class_distribution.gex.svg')
        pd.DataFrame(columns=["class", "Count"]).to_csv('classes.gex.tsv', sep='\t', index=False)

    # Handle DNAM
    if has_dnam.any():
        handle_modality(has_dnam, 'DNAM', 'classes.dnam.tsv', 'class_distribution.dnam.svg')
    else:
        plot_empty('No DNAM samples found', 'Class Distribution (DNAM)', 'class_distribution.dnam.svg')
        pd.DataFrame(columns=["class", "Count"]).to_csv('classes.dnam.tsv', sep='\t', index=False)

    # Create versions file
    import numpy
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pd.__version__}"\\n')
        f.write(f'    numpy: "{numpy.__version__}"\\n')
        f.write(f'    matplotlib: "{matplotlib.__version__}"\\n')
    """

    stub:
    """
    touch class_distribution.gex.svg
    touch class_distribution.dnam.svg
    touch classes.gex.tsv
    touch classes.dnam.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
