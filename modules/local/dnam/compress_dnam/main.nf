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

if "probe_id" in dataframe.columns:
    dataframe = dataframe.set_index("probe_id")
else:
    dataframe = dataframe.set_index(dataframe.columns[0])

dataframe = dataframe.apply(pd.to_numeric, errors="coerce").round(4).astype("float32")
dataframe.to_csv("compressed_betas.tsv", sep="\t", index_label="probe_id")

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
