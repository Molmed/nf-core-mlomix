process CLASS_REPORTER {
    tag "class_reporter"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/pandas:1.5.3_cv1"

    input:
    path samplesheet

    output:
    path "class_distribution.svg"  , emit: class_plot
    path "class_distribution.csv"  , emit: class_csv
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

    def plot_class_distribution(y_series, title):
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

    def plot_empty(message, title):
        plt.figure(figsize=(10, 5))
        plt.text(0.5, 0.5, message, ha='center', va='center', fontsize=14, transform=plt.gca().transAxes)
        plt.title(title)
        plt.xlabel("Class")
        plt.ylabel("Count")
        plt.savefig("class_distribution.svg", bbox_inches='tight')

    # Auto-detect delimiter
    samplesheet_path = "${samplesheet}"
    if samplesheet_path.endswith('.tsv'):
        sep = '\t'
    else:
        sep = ','

    df = pd.read_csv(samplesheet_path, sep=sep)

    if 'class' not in df.columns:
        print("WARNING: Samplesheet does not contain a 'class' column. Generating empty report.")
        plot_empty("No 'class' column found in samplesheet", "Class Distribution")
        pd.DataFrame(columns=["class", "Count"]).to_csv("class_distribution.csv", index=False)
    else:
        classes = df['class'].dropna()
        if classes.empty:
            print("WARNING: 'class' column is present but contains no values. Generating empty report.")
            plot_empty("'class' column has no values", "Class Distribution")
            pd.DataFrame(columns=["class", "Count"]).to_csv("class_distribution.csv", index=False)
        else:
            plot_class_distribution(classes, "Class Distribution")
            classes.value_counts().sort_values(ascending=False).to_csv(
                "class_distribution.csv", header=["Count"]
            )

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
    touch class_distribution.svg
    touch class_distribution.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
