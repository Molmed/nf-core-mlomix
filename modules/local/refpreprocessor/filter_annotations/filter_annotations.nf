process FILTER_ANNOTATIONS {
    tag "$annotations"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path annotations
    val full_genome_name

    output:
    path "${full_genome_name}.annotations.filtered.csv", emit: filtered_annotations
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_file = "${full_genome_name}.annotations.filtered.csv"

    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    def filter_annotations(input_file, output_file):
        # Read in file
        annot = pd.read_csv(input_file)

        # Get all rows where biotype is protein_coding
        protein_coding_annot = annot.loc[annot['biotype'] == 'protein_coding']

        # Only keep genes that are on chrs 1-22 and X
        valid_chromosomes = [str(i) for i in range(1, 23)] + ['X']

        # Keep only genes where chr is in valid_chromosomes
        filtered_annot = protein_coding_annot[
            protein_coding_annot['chr'].isin(valid_chromosomes)]

        # Replace NA in names with empty string
        filtered_annot['name'] = filtered_annot['name'].fillna('')

        # Remove ribosomal genes
        filtered_annot = filtered_annot[
            ~filtered_annot['name'].str.startswith(('RPL', 'RPS'))]

        # Write the filtered DataFrame to a CSV file
        filtered_annot.to_csv(output_file, index=False)


    # Filter the annotations
    filter_annotations("${annotations}", "${output_file}")

    # Create versions file
    import pandas
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
    """
}
