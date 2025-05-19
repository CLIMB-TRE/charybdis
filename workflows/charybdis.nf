/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_charybdis_pipeline'

include { ONT_ASSEMBLY           } from '../subworkflows/local/ont_assembly/main'
include { ILLUMINA_ASSEMBLY      } from '../subworkflows/local/illumina_assembly/main'

include { KRAKEN2_KRAKEN2        } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_CLIENT         } from '../modules/local/kraken2-client/main'
include { METABAT2_METABAT2      } from '../modules/nf-core/metabat2/metabat2/main'
include { BANDAGE_IMAGE          } from '../modules/nf-core/bandage/image/main'
include { UNTAR                  } from '../modules/nf-core/untar/main'
include { AMRFINDERPLUS_UPDATE   } from '../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN      } from '../modules/nf-core/amrfinderplus/run/main'

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
            illumina: meta.platform == "illumina" || meta.platform == "illumina.se"
        }
        .set { ch_input }

    ONT_ASSEMBLY(ch_input.ont)
    ch_versions = ch_versions.mix(ONT_ASSEMBLY.out.versions.first())

    ILLUMINA_ASSEMBLY(ch_input.illumina)
    ch_versions = ch_versions.mix(ILLUMINA_ASSEMBLY.out.versions.first())

    ch_contigs = ONT_ASSEMBLY.out.contigs.mix(ILLUMINA_ASSEMBLY.out.contigs)
    ch_graph = ONT_ASSEMBLY.out.gfa.mix(ILLUMINA_ASSEMBLY.out.fastg)

    if (!params.k2_remote) {

        if (!params.k2_local) {
            error("No Kraken2 database provided. Please provide a local or remote Kraken2 database.")
        }

        // Check if the local Kraken2 database is a tarball
        // If it is, extract it and set the path to the extracted database
        if (params.k2_local.endsWith(".tar.gz") || params.k2_local.endsWith(".tgz")) {
            ch_k2_db_untar = file(params.k2_local)
            UNTAR([[:], ch_k2_db_untar])
            ch_versions = ch_versions.mix(UNTAR.out.versions)

            UNTAR.out.untar
                .map { _meta, path -> path }
                .set { ch_k2_db_local }
        }
        else {
            ch_k2_db_local = file(params.k2_local)
        }

        // Even if input was paired, contigs are always single-end
        ch_contigs.map { meta, contigs -> [[id: meta.id, single_end: true], contigs] }.set { ch_k2_local_input }
        KRAKEN2_KRAKEN2(
            ch_k2_local_input,
            ch_k2_db_local,
            false,
            true,
        )
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())
    }
    else {
        KRAKEN2_CLIENT(
            ch_contigs,
            params.k2_remote,
        )
        ch_versions = ch_versions.mix(KRAKEN2_CLIENT.out.versions.first())
    }

    // Run AMRFinderPlus
    AMRFINDERPLUS_UPDATE()
    ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)

    AMRFINDERPLUS_RUN(ch_contigs, AMRFINDERPLUS_UPDATE.out.db)
    ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions.first())

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
