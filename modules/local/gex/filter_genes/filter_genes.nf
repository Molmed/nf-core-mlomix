process FILTER_GENES {
    label 'process_single'
    scratch true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_gene-thesaurus:587dd866009f8a59' :
        'community.wave.seqera.io/library/pip_gene-thesaurus:7c78bdcbd4b116b6' }"

    input:
    val dataset_name
    path data_path
    each path(ref_path)

    output:
    path "${dataset_name}.filtered_genes.csv", emit: filtered_genes_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gt_cache_dir = "${workDir}/gene_thesaurus_cache"
    """
    #!/usr/bin/env python3

    import os
    import pandas as pd
    from gene_thesaurus import GeneThesaurus

    def filter():
        data = pd.read_csv("${data_path}", index_col=0)
        ref = pd.read_csv("${ref_path}")

        # Use shared cache directory in Nextflow work dir
        cache_dir = "${gt_cache_dir}"
        os.makedirs(cache_dir, exist_ok=True)
        gt = GeneThesaurus(data_dir=cache_dir)

        # Save all the sample columns
        case_columns = data.columns

        # Determine input gene format
        # Take the first gene. If it begins with "ENS", it is "ensembl"
        # If it is an integer, it is "entrez"
        # Otherwise it is "symbol"

        if data.index[0].startswith("ENS"):
            gene_format = 'ensembl_id'
        elif data.index[0].isdigit():
            gene_format = 'entrez_id'
        else:
            gene_format = 'symbol'

        if gene_format != 'symbol':
            # Get gene names from ensembl ids
            translated_genes = gt.translate_genes(data.index.values,
                                                        source=gene_format,
                                                        target='symbol')

            def _standardize_gene_name(gene_name):
                if gene_name in translated_genes:
                    return translated_genes[gene_name]
                return gene_name

            data['gene_name_std'] = data.index.map(_standardize_gene_name)
        else:
            # Code to use if index is in symbol format
            # Get all gene name values to update to the latest standard
            updated_genes = gt.update_gene_symbols(data.index.values)

            # Copy data.index to gene_name_std
            data['gene_name_std'] = data.index

            # Update gene_name_std from updated_genes
            data['gene_name_std'] = data['gene_name_std'].map(
                updated_genes).fillna(data['gene_name_std'])

        # Do the same for ref
        updated_genes = gt.update_gene_symbols(ref['name'].values)

        # Set new column with the latest gene name
        ref['gene_name_std'] = ref['name'].map(
            updated_genes).fillna(ref['name'])

        filtered_rows = {}

        def _find_ref_key(row):
            if row.name in ref['id'].values:
                return ref[ref['id'] == row.name]['id'].values[0]
            if row['gene_name_std'] in ref['gene_name_std'].values:
                ref_rows = ref[ref['gene_name_std'] == row['gene_name_std']]

                if len(ref_rows) == 1:
                    return ref_rows['id'].values[0]

                # See if one of the rows also matches a ref id
                for i, r in ref_rows.iterrows():
                    if r['id'] in ref['id'].values:
                        return r['id']

                # Return the first row
                return ref_rows['id'].values[0]

        for i, row in data.iterrows():
            # Strip whitespace
            key = _find_ref_key(row)

            # Print ref row length
            if key is not None:
                try:
                    key = key.strip()
                except Exception:
                    print(f'Error: {key}')
                    exit()

                # Just keep case cols in row
                row = row[case_columns]

                # Is there already a row in filtered_rows?
                if key in filtered_rows:
                    filtered_rows[key] += row
                else:
                    filtered_rows[key] = row
                # Keep the previous index name for the row
                filtered_rows[key].name = key

        data = pd.DataFrame(filtered_rows.values())

        # Print records in ref.id that are not in data.id
        missing = ref[~ref['id'].isin(data.index)]

        # Create records for all missing genes in data, filled with 0s
        # The missing id is the index value, and all the case columns are 0
        missing_data = pd.DataFrame(index=missing['id'],
                                    columns=case_columns,
                                    data=0)

        # Remove index name
        missing_data.index.name = None

        # Append the missing data to the data
        data = pd.concat([data, missing_data])

        # Sort data by index
        data = data.sort_index()

        # Write to file
        data.to_csv("${dataset_name}.filtered_genes.csv")

    # Filter the annotations
    filter()

    # Create versions file
    import pandas
    import sys
    from importlib.metadata import version

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
        f.write(f'    gene-thesaurus: "{version("gene-thesaurus")}"\\n')
    """
}
