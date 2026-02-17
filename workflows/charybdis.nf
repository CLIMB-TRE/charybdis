/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap              } from 'plugin/nf-schema'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../subworkflows/local/utils_nfcore_charybdis_pipeline'

include { ONT_ASSEMBLY                  } from '../subworkflows/local/ont_assembly/main'
include { ILLUMINA_ASSEMBLY             } from '../subworkflows/local/illumina_assembly/main'

include { KRAKEN2_KRAKEN2               } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_CLIENT                } from '../modules/local/kraken2-client/main'
include { METABAT2_METABAT2             } from '../modules/nf-core/metabat2/metabat2/main'
include { GTDBTK_CLASSIFYWF             } from '../modules/nf-core/gtdbtk/classifywf/main'
include { GTDBTK_GTDBTONCBIMAJORITYVOTE } from '../modules/nf-core/gtdbtk/gtdbtoncbimajorityvote/main'
include { BANDAGE_IMAGE                 } from '../modules/nf-core/bandage/image/main'
include { UNTAR as UNTAR_KRAKEN         } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_TAXONOMY       } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_GTDB           } from '../modules/nf-core/untar/main'
include { AMRFINDERPLUS_UPDATE          } from '../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN             } from '../modules/nf-core/amrfinderplus/run/main'
include { TAXONKIT_LINEAGE              } from '../modules/nf-core/taxonkit/lineage/main'
include { MINIMAP2_ALIGN                } from '../modules/nf-core/minimap2/align/main'
include { BWAMEM2_INDEX                 } from '../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM                   } from '../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX                } from '../modules/nf-core/samtools/index/main'
include { DIVERSITY_METRICS             } from '../modules/local/diversity-metrics/main'
include { PERBASE                       } from '../modules/nf-core/perbase/main'
include { CALCULATE_PER_TAXON_DIVERSITY } from '../modules/local/per-taxon-diversity/main'
include { DIVERSITY_METRIC_PLOT         } from '../modules/local/diversity-metric-plot/main'
include { BUSCO_BUSCO                   } from '../modules/nf-core/busco/busco/main'
include { GENOMAD_DOWNLOAD              } from '../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND              } from '../modules/nf-core/genomad/endtoend/main'

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
    ch_graph = ONT_ASSEMBLY.out.gfa.mix(ILLUMINA_ASSEMBLY.out.graph)

    if (!params.k2_remote) {

        if (!params.k2_local) {
            error("No Kraken2 database provided. Please provide a local or remote Kraken2 database.")
        }

        // Check if the local Kraken2 database is a tarball
        // If it is, extract it and set the path to the extracted database
        if (params.k2_local.endsWith(".tar.gz") || params.k2_local.endsWith(".tgz")) {
            k2_db_tarball = file(params.k2_local, checkIfExists: true)
            UNTAR_KRAKEN([[:], k2_db_tarball])
            ch_versions = ch_versions.mix(UNTAR_KRAKEN.out.versions)

            UNTAR_KRAKEN.out.untar
                .map { _meta, path -> path }
                .set { k2_db_local }
        }
        else {
            k2_db_local = file(params.k2_local, checkIfExists: true)
        }

        // Even if input was paired, contigs are always single-end -> don't modify
        ch_k2_local_input = ch_contigs.map { meta, contigs -> [(meta - meta.subMap("single_end")) + [single_end: true], contigs] }


        KRAKEN2_KRAKEN2(
            ch_k2_local_input,
            k2_db_local,
            false,
            true,
        )
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

        ch_k2_report = KRAKEN2_KRAKEN2.out.report
        ch_k2_assignments = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    }
    else {
        KRAKEN2_CLIENT(
            ch_contigs,
            params.k2_remote,
        )
        ch_versions = ch_versions.mix(KRAKEN2_CLIENT.out.versions.first())

        ch_k2_report = KRAKEN2_CLIENT.out.report
        ch_k2_assignments = KRAKEN2_CLIENT.out.assignments
    }

    // If a taxonomy file is provided as a tarball, extract it
    if (params.taxonomy.endsWith(".tar.gz") || params.taxonomy.endsWith(".tgz")) {
        taxonomy_tarball = file(params.taxonomy, checkIfExists: true)
        UNTAR_TAXONOMY([[:], taxonomy_tarball])
        ch_versions = ch_versions.mix(UNTAR_TAXONOMY.out.versions)

        UNTAR_TAXONOMY.out.untar
            .map { _meta, path -> path }
            .set { taxonomy }
    }
    else {
        taxonomy = file(params.taxonomy, checkIfExists: true)
    }

    // Run taxonkit lineage on the Kraken2 assignments
    TAXONKIT_LINEAGE(
        ch_k2_assignments.map { meta, assignments -> [meta, [], assignments] },
        taxonomy,
    )
    ch_versions = ch_versions.mix(TAXONKIT_LINEAGE.out.versions.first())

    ch_contigs_branched = ch_contigs.branch { meta, _contigs ->
        ont: meta.platform == "ont"
        illumina: meta.platform == "illumina" || meta.platform == "illumina.se"
    }

    BWAMEM2_INDEX(ch_contigs_branched.illumina)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

    ch_bwamem2_input = ch_input.illumina
        .join(BWAMEM2_INDEX.out.index)
        .join(ch_contigs_branched.illumina)
        .multiMap { meta, fastq_1, fastq_2, contig_index, contigs ->
            reads: [meta, [fastq_1, fastq_2]]
            contig_index: [meta, contig_index]
            contig_fasta: [meta, contigs]
        }

    BWAMEM2_MEM(
        ch_bwamem2_input.reads,
        ch_bwamem2_input.contig_index,
        ch_bwamem2_input.contig_fasta,
        true,
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    SAMTOOLS_INDEX(BWAMEM2_MEM.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_minimap2_input = ch_input.ont
        .join(ch_contigs_branched.ont)
        .multiMap { meta, reads1, _reads2, contigs ->
            reads: [meta, reads1]
            contig_fasta: [meta, contigs]
        }

    MINIMAP2_ALIGN(ch_minimap2_input.reads, ch_minimap2_input.contig_fasta, true, "bai", true, false)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    ch_bams_with_index = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index)
        .mix(
            BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)
        )
        .map { meta, bam, bai ->
            [meta, bam, bai, []]
        }

    PERBASE(ch_bams_with_index, [[:], [], []])
    ch_versions = ch_versions.mix(PERBASE.out.versions.first())

    // Generate diversity metrics for the contigs
    // These processes literally just use Python standard library, so no need for versions files (nobody cares about Python versions)
    DIVERSITY_METRICS(PERBASE.out.tsv)

    ch_lineages_and_diversity = TAXONKIT_LINEAGE.out.tsv
        .map { meta, tsv -> [meta.subMap("id", "platform"), tsv] }
        .join(
            DIVERSITY_METRICS.out.csv.map { meta, csv -> [meta.subMap("id", "platform"), csv] }
        )

    CALCULATE_PER_TAXON_DIVERSITY(ch_lineages_and_diversity)

    DIVERSITY_METRIC_PLOT(
        CALCULATE_PER_TAXON_DIVERSITY.out.taxon_diversity_tsv
    )

    // Run AMRFinderPlus
    AMRFINDERPLUS_UPDATE()
    ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)

    AMRFINDERPLUS_RUN(ch_contigs, AMRFINDERPLUS_UPDATE.out.db)
    ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)

    // Generate a Bandage image of the assembly graph (it says it requires GFA but works fine with fastg)
    BANDAGE_IMAGE(
        ch_graph
    )
    ch_versions = ch_versions.mix(BANDAGE_IMAGE.out.versions.first())

    // Bin the contigs with metabat2
    if (!params.skip_binning) {
        METABAT2_METABAT2(
            ch_contigs.map { meta, contigs -> [meta, contigs, []] }
        )
        ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())
        ch_bins = METABAT2_METABAT2.out.fasta
    }
    else {
        ch_bins = Channel.empty()
    }


    //
    // Prep the GTDB database
    //
    if (params.run_gtdb_tk) {
        if (params.gtdb_db.endsWith(".tar.gz") || params.gtdb_db.endsWith(".tgz")) {
            gtdb_tarball = file(params.gtdb_db, checkIfExists: true)
            UNTAR_GTDB([[:], gtdb_tarball])
            ch_versions = ch_versions.mix(UNTAR_GTDB.out.versions)

            gtdb_db = UNTAR_GTDB.out.untar.map { _meta, path -> [params.gtdb_version, path] }
        }
        else {
            gtdb_val = file(params.gtdb_db, checkIfExists: true)
            gtdb_db = [gtdb_val.simpleName, gtdb_val]
        }

        GTDBTK_CLASSIFYWF(
            ch_bins,
            gtdb_db,
            false,
        )
        ch_versions = ch_versions.mix(GTDBTK_CLASSIFYWF.out.versions.first())

    }


    // GTDBTK_GTDBTONCBIMAJORITYVOTE(
    //     ch_gtdb_convert_input,
    //     gtdb_ar53,
    //     gtdb_bac120,
    // )
    // ch_versions = ch_versions.mix(GTDBTK_GTDBTONCBIMAJORITYVOTE.out.versions.first())

    //
    // Assess bin quality with BUSCO
    //
    store_dir = file(params.store_dir, checkIfExists: true)
    busco_store_dir = file(
        "${store_dir.toUriString()}/busco_db"
    )

    // Try to figure out how to store the busco db between runs later
    BUSCO_BUSCO(
        ch_bins,
        "genome",
        "auto",
        [],
        [],
        true,
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // Search for mobile elements (plasmids, viruses, etc) with genomad
    //
    GENOMAD_DOWNLOAD()
    ch_versions = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions.first())

    GENOMAD_ENDTOEND(
        ch_contigs,
        GENOMAD_DOWNLOAD.out.genomad_db,
    )
    ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions.first())

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
