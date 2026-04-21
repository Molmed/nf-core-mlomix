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

workflow GEX_REF_PREPROCESSOR {

    take:
    genome             //  string: Version of the genome to use, e.g. 'GRCh38' or 'hg38'
    annotation_version //  int: Version of the annotation to use, e.g. '109' for Ensembl release 109

    main:

    ch_versions = channel.empty()
    ch_missing_genome = channel.empty()
    ch_missing_annotation_version = channel.empty()
    ch_full_genome_name_missing = channel.empty()

    ch_cached_filtered_annotations = genome
        .combine(annotation_version)
        .map { g, ann -> file("${params.annotation_cache_dir}/Homo_sapiens.${g}.${ann}.annotations.filtered.csv") }
        .filter { cached -> cached.exists() }

    ch_missing_refs = genome
        .combine(annotation_version)
        .filter { g, ann -> !file("${params.annotation_cache_dir}/Homo_sapiens.${g}.${ann}.annotations.filtered.csv").exists() }

    ch_missing_genome = ch_missing_refs.map { g, _ann -> g }
    ch_missing_annotation_version = ch_missing_refs.map { _g, ann -> ann }
    ch_full_genome_name_missing = ch_missing_refs.map { g, ann -> "Homo_sapiens.${g}.${ann}" }

    //
    // MODULE: Download GTF file
    //
    DOWNLOAD_GTF (
        ch_missing_genome,
        ch_missing_annotation_version
    )
    ch_versions = ch_versions.mix(DOWNLOAD_GTF.out.versions)

    //
    // MODULE: Flatten GTF file
    //
    // TODO: Rename gtf_file to gtf
    FLATTEN_GTF (
        DOWNLOAD_GTF.out.gtf_file,
        ch_full_genome_name_missing
    )
    ch_versions = ch_versions.mix(FLATTEN_GTF.out.versions.first())

    //
    // MODULE: Parse GTF file
    //
    PARSE_GTF (
        DOWNLOAD_GTF.out.gtf_file,
        FLATTEN_GTF.out.saf,
        ch_full_genome_name_missing
    )
    ch_versions = ch_versions.mix(PARSE_GTF.out.versions.first())

    //
    // MODULE: Filter annotations
    //
    FILTER_ANNOTATIONS (
        PARSE_GTF.out.annotations,
        ch_full_genome_name_missing
    )
    ch_versions = ch_versions.mix(FILTER_ANNOTATIONS.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'gex_ref_preprocessor_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    filtered_annotations = ch_cached_filtered_annotations.mix(FILTER_ANNOTATIONS.out.filtered_annotations) // channel: [ path(annotations.filtered.csv) ]
    versions           = ch_versions                          // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
