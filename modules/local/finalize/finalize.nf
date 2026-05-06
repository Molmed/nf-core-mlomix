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
    path visual_files

    output:
    path "labels.csv", emit: labels_csv
    path "features.dnam.csv", emit: dnam_features_csv
    path "features.gex.csv", emit: gex_features_csv
    path "feature_names.dnam.txt", emit: dnam_feature_names
    path "feature_names.gex.txt", emit: gex_feature_names
    path "visuals", emit: visuals_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def visual_file_names = (visual_files instanceof List ? visual_files : [visual_files])
        .collect { visual_file -> "\"${visual_file.name}\"" }
        .join(', ')
"""
#!/usr/bin/env python3

import pandas as pd
import sys
import os
import shutil

# Read labels (classes)
labels = pd.read_csv("${classes_file}", sep="\t", index_col=0)
labels = labels.iloc[:, 0]
label_samples = labels.index.tolist()

# Write labels as labels.csv with index renamed to "id"
labels.index.name = "id"
labels.to_csv("labels.csv", sep=",")

# Read and reindex DNAM features (only if file provided and not empty)
try:
    dnam = pd.read_csv("${dnam_transposed}", sep=",", index_col=0)
    dnam = dnam.loc[label_samples]
    if not dnam.empty:
        dnam.to_csv("features.dnam.csv", sep=",")
    else:
        pd.DataFrame().to_csv("features.dnam.csv", sep=",")
except FileNotFoundError:
    pd.DataFrame().to_csv("features.dnam.csv", sep=",")  # Create empty file

# Read and reindex GEX features (only if file provided and not empty)
try:
    gex = pd.read_csv("${gex_transposed}", sep=",", index_col=0)
    gex = gex.loc[label_samples]
    if not gex.empty:
        gex.to_csv("features.gex.csv", sep=",")
    else:
        pd.DataFrame().to_csv("features.gex.csv", sep=",")
except FileNotFoundError:
    pd.DataFrame().to_csv("features.gex.csv", sep=",")  # Create empty file

# Write feature names for DNAM
try:
    dnam = pd.read_csv("${dnam_transposed}", sep=",", index_col=0)
    dnam = dnam.loc[label_samples]
    if not dnam.empty:
        with open("feature_names.dnam.txt", "w") as f:
            f.write("\\n".join(dnam.columns.tolist()))
    else:
        open("feature_names.dnam.txt", "w").close()
except FileNotFoundError:
    open("feature_names.dnam.txt", "w").close()

# Write feature names for GEX
try:
    gex = pd.read_csv("${gex_transposed}", sep=",", index_col=0)
    gex = gex.loc[label_samples]
    if not gex.empty:
        with open("feature_names.gex.txt", "w") as f:
            f.write("\\n".join(gex.columns.tolist()))
    else:
        open("feature_names.gex.txt", "w").close()
except FileNotFoundError:
    open("feature_names.gex.txt", "w").close()
os.makedirs("visuals", exist_ok=True)
visual_inputs = [${visual_file_names}]
for visual_path in visual_inputs:
    if not visual_path.lower().endswith((".png", ".svg")):
        continue

    destination = os.path.join("visuals", os.path.basename(visual_path))
    if not os.path.exists(destination):
        shutil.copy2(visual_path, destination)
        continue

    stem, ext = os.path.splitext(os.path.basename(visual_path))
    idx = 1
    while True:
        candidate = os.path.join("visuals", f"{stem}.{idx}{ext}")
        if not os.path.exists(candidate):
            shutil.copy2(visual_path, candidate)
            break
        idx += 1

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
    touch feature_names.dnam.txt
    touch feature_names.gex.txt
    mkdir -p visuals

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        pandas: "stub"
    END_VERSIONS
    """
}
