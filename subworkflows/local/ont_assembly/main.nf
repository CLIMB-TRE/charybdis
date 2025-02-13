include { FLYE              } from '../../../modules/nf-core/flye/main'
include { PORECHOP_PORECHOP } from '../../../modules/nf-core/porechop/porechop/main'
include { PIGZ_UNCOMPRESS   } from '../../../modules/nf-core/pigz/uncompress/main'

workflow ONT_ASSEMBLY {
    take:
    ch_input

    main:

    ch_versions = Channel.empty()

    PORECHOP_PORECHOP(
        ch_input.map { meta, fastq_1, _fastq_2 -> [meta, fastq_1] }
    )
    ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())

    // Add some logic to determine the mode?
    mode = "--nano-corr"
    FLYE(
        PORECHOP_PORECHOP.out.reads,
        mode,
    )
    ch_versions = ch_versions.mix(FLYE.out.versions.first())

    PIGZ_UNCOMPRESS(
        FLYE.out.gfa
    )
    ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS.out.versions.first())

    emit:
    reads    = PORECHOP_PORECHOP.out.reads
    contigs  = FLYE.out.fasta
    gfa      = PIGZ_UNCOMPRESS.out.file
    versions = ch_versions // channel: [ versions.yml ]
}
