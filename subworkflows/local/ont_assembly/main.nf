include { FLYE              } from '../../../modules/nf-core/flye/main'
include { PORECHOP_PORECHOP } from '../../../modules/nf-core/porechop/porechop/main'

workflow ONT_ASSEMBLY {
    take:
    ch_input

    main:

    ch_versions = Channel.empty()

    ch_fastq = ch_input.map { meta, platform, fastq_1, _fastq_2 -> [[id: meta.id, platform: platform, single_end: true], fastq_1] }

    PORECHOP_PORECHOP(
        ch_fastq
    )
    ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())

    // Add some logic to determine the mode?
    mode = "--nano-corr"
    FLYE(
        PORECHOP_PORECHOP.out.reads,
        mode,
    )
    ch_versions = ch_versions.mix(FLYE.out.versions.first())

    emit:
    reads    = PORECHOP_PORECHOP.out.reads
    contigs  = FLYE.out.fasta
    gfa      = FLYE.out.gfa
    versions = ch_versions // channel: [ versions.yml ]
}
