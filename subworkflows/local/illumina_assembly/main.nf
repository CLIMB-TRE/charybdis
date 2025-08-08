include { MEGAHIT       } from '../../../modules/nf-core/megahit/main'
include { MEGAHIT_FASTG } from '../../../modules/local/megahit/fastg/main'
include { SPADES        } from '../../../modules/nf-core/spades/main'

workflow ILLUMINA_ASSEMBLY {
    take:
    ch_input

    main:

    ch_versions = Channel.empty()

    if (params.illumina_assembler == "spades") {
        ch_spades_input = ch_input.map { meta, fastq1, fastq2 ->
            [meta, [fastq1, fastq2], [], []]
        }

        SPADES(ch_spades_input, [], [])
        ch_versions = ch_versions.mix(SPADES.out.versions.first())

        ch_contigs = SPADES.out.contigs
        ch_graph = SPADES.out.gfa
    }
    else if (params.illumina_assembler == "megahit") {
        MEGAHIT(ch_input)
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())

        MEGAHIT_FASTG(MEGAHIT.out.k_contigs)
        ch_versions = ch_versions.mix(MEGAHIT_FASTG.out.versions.first())

        ch_contigs = MEGAHIT.out.contigs
        ch_graph = MEGAHIT_FASTG.out.fastg
    }
    else {
        error("Unknown Illumina assembler: ${params.illumina_assembler}")
    }

    emit:
    contigs  = ch_contigs
    graph    = ch_graph
    versions = ch_versions // channel: [ path(versions.yml) ]
}
