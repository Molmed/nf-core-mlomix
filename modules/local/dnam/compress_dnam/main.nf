process COMPRESS_DNAM {
    tag "compress_dnam"
    label 'process_medium'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(dataset_name), val(sample_name), path(betas)

    output:
    tuple val(dataset_name), val(sample_name), path("${dataset_name}__${sample_name}.compressed_betas.tsv"), emit: compressed_betas
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_name = "${dataset_name}__${sample_name}.compressed_betas.tsv"
    """
    python3 <<'PY'
import pandas as pd
import sys

dataframe = pd.read_csv("${betas}", sep="\t")

if dataframe.shape[1] < 2:
    raise ValueError(
        f"Expected at least 2 columns in ${betas}, got {dataframe.shape[1]}"
    )

first_col = str(dataframe.columns[0])
looks_like_probe_id = first_col.lower().startswith(("cg", "ch", "rs"))

if first_col != "probe_id" and looks_like_probe_id:
    dataframe = pd.read_csv("${betas}", sep="\t", header=None)

if str(dataframe.iloc[0, 0]).strip().lower() == "probe_id":
    dataframe = dataframe.iloc[1:, :]

dataframe.columns = ["probe_id"] + [f"value_{i}" for i in range(1, dataframe.shape[1])]
dataframe = dataframe.set_index("probe_id")

dataframe = dataframe.apply(pd.to_numeric, errors="coerce").round(4).astype("float32")
dataframe.index.name = None
dataframe.to_csv("compressed_betas.tsv", sep="\t", header=False)

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'    python: "{sys.version.split()[0]}"\\n')
    f.write(f'    pandas: "{pd.__version__}"\\n')
PY

    mv compressed_betas.tsv ${output_name}
    """

    stub:
    """
    touch ${dataset_name}__${sample_name}.compressed_betas.tsv
    printf '"%s":\n    python: "stub"\n    pandas: "stub"\n' "${task.process}" > versions.yml
    """
}
