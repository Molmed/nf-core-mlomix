/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_GTF           } from '../modules/local/refpreprocessor/download_gtf'
include { FLATTEN_GTF            } from '../modules/local/refpreprocessor/flatten_gtf'
include { PARSE_GTF              } from '../modules/local/refpreprocessor/parse_gtf'
include { FILTER_ANNOTATIONS     } from '../modules/local/refpreprocessor/filter_annotations/filter_annotations'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REFPREPROCESSOR {

    take:
    genome             //  string: Version of the genome to use, e.g. 'GRCh38' or 'hg38'
    annotation_version //  int: Version of the annotation to use, e.g. '109' for Ensembl release 109

    main:

    ch_versions = channel.empty()
    full_genome_name = "Homo_sapiens.${genome}.${annotation_version}"

    //
    // MODULE: Download GTF file
    //
    DOWNLOAD_GTF (
        genome,
        annotation_version
    )
    ch_versions = ch_versions.mix(DOWNLOAD_GTF.out.versions)

    //
    // MODULE: Flatten GTF file
    //
    // TODO: Rename gtf_file to gtf
    FLATTEN_GTF (
        DOWNLOAD_GTF.out.gtf_file,
        full_genome_name
    )
    ch_versions = ch_versions.mix(FLATTEN_GTF.out.versions.first())

    //
    // MODULE: Parse GTF file
    //
    PARSE_GTF (
        DOWNLOAD_GTF.out.gtf_file,
        FLATTEN_GTF.out.saf,
        full_genome_name
    )
    ch_versions = ch_versions.mix(PARSE_GTF.out.versions.first())

    //
    // MODULE: Filter annotations
    //
    FILTER_ANNOTATIONS (
        PARSE_GTF.out.annotations,
        full_genome_name
    )
    ch_versions = ch_versions.mix(FILTER_ANNOTATIONS.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'refpreprocessor_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    filtered_annotations = FILTER_ANNOTATIONS.out.filtered_annotations // channel: [ path(annotations.filtered.csv) ]
    versions           = ch_versions                          // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
