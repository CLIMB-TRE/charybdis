/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                              } from 'plugin/nf-schema'
include { softwareVersionsToYAML                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                        } from '../subworkflows/local/utils_nfcore_charybdis_pipeline'

include { ONT_ASSEMBLY                                  } from '../subworkflows/local/ont_assembly/main'
include { ILLUMINA_ASSEMBLY as ILLUMINA_ASSEMBLY_PAIRED } from '../subworkflows/local/illumina_assembly/main'
include { ILLUMINA_ASSEMBLY as ILLUMINA_ASSEMBLY_SINGLE } from '../subworkflows/local/illumina_assembly/main'

include { KRAKEN2_KRAKEN2                               } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_CLIENT                                } from '../modules/local/kraken2-client/main'
include { METABAT2_METABAT2                             } from '../modules/nf-core/metabat2/metabat2/main'
include { BANDAGE_IMAGE                                 } from '../modules/nf-core/bandage/image/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHARYBDIS {
    take:
    ch_samplesheet

    main:

    ch_versions = Channel.empty()

    // Run the appropriate assembly pipeline for the platform
    ch_samplesheet
        .branch { meta, _fastq_1, _fastq_2 ->
            ont: meta.platform == "ont"
            illumina: meta.platform == "illumina"
            illumina_se: meta.platform == "illumina.se"
        }
        .set { ch_input }

    ONT_ASSEMBLY(ch_input.ont)
    ch_versions = ch_versions.mix(ONT_ASSEMBLY.out.versions.first())

    ILLUMINA_ASSEMBLY_PAIRED(ch_input.illumina)
    ch_versions = ch_versions.mix(ILLUMINA_ASSEMBLY_PAIRED.out.versions.first())

    ILLUMINA_ASSEMBLY_SINGLE(ch_input.illumina_se)
    ch_versions = ch_versions.mix(ILLUMINA_ASSEMBLY_SINGLE.out.versions.first())

    ch_contigs = ONT_ASSEMBLY.out.contigs.mix(ILLUMINA_ASSEMBLY_PAIRED.out.contigs, ILLUMINA_ASSEMBLY_SINGLE.out.contigs)
    ch_graph = ONT_ASSEMBLY.out.gfa.mix(ILLUMINA_ASSEMBLY_PAIRED.out.fastg, ILLUMINA_ASSEMBLY_SINGLE.out.fastg)

    if (!params.k2_remote) {
        ch_contigs.map { meta, contigs -> [[id: meta.id, single_end: true], contigs] }.set { ch_k2_local_input }
        KRAKEN2_KRAKEN2(
            ch_k2_local_input,
            params.k2_local,
            false,
            true,
        )
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())
    }
    else {
        // Even if input was paired, contigs are always single-end
        KRAKEN2_CLIENT(
            ch_contigs,
            params.k2_remote,
        )
        ch_versions = ch_versions.mix(KRAKEN2_CLIENT.out.versions.first())
    }

    // Generate a Bandage image of the assembly graph (it says it requires GFA but works fine with fastg)
    BANDAGE_IMAGE(
        ch_graph
    )
    ch_versions = ch_versions.mix(BANDAGE_IMAGE.out.versions.first())

    // Bin the contigs with metabat2
    METABAT2_METABAT2(
        ch_contigs.map { meta, contigs -> [meta, contigs, []] }
    )
    ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'charybdis_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }
}
