# nf-core/mlomix pipeline - params.input schema

Schema for the file provided with params.input

| Column | Description | Type | Required | Pattern |
|---|---|---|---|---|
| `sample` | Unique sample identifier, without spaces. | string | Yes | `^\S+$` |
| `dataset` | Unique dataset identifier, without spaces, e.g. sequencing run id, project identifier, etc. | string | Yes | `^\S+$` |
| `batch` | Batch identifier, without spaces. Only required for model training with samples from multiple sequencing batches. | string | No | `^\S+$` |
| `class` | Class/subtype label for the sample. Only required for model training or validation. | string | No |  |
| `data_silo` | Data silo identifier, e.g. 'train', 'test', 'validation'. Only used for model training. | string | No |  |
| `gex_feature_counts_file` | Path to the gene expression feature counts file, if available for this sample. The file should be tab-delimited with the first column containing gene identifiers and the second column containing raw counts. | file-path | No |  |
| `dnam_beta_matrix_file` | Path to the DNA methylation beta matrix file, if available for this sample. The file should be tab-delimited with the first column containing probe identifiers and the second column containing beta values. | file-path | No |  |
| `dnam_pvals_file` | Path to the DNA methylation p-values file, if available for this sample. The file should be tab-delimited with the first column containing probe identifiers and the second column containing p-values. | file-path | No |  |
| `sentrix_id` | Sentrix ID, e.g. the ID of the array slide. Only required if processing raw DNA methylation data from Illumina arrays (IDAT). | ['string', 'integer'] | No |  |
| `sentrix_position` | Sentrix position, e.g. the position of the array slide. Only required if processing raw DNA methylation data from Illumina arrays (IDAT). | string | No | `^\S+$` |
| `idats_basename` | Path to directory containing *_Red.idat and *_Grn.idat files. Only required if processing raw DNA methylation data from Illumina arrays (IDAT). | directory-path | No |  |
