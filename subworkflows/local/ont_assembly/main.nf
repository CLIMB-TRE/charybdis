include { FLYE              } from '../../../modules/nf-core/flye/main'
include { PORECHOP_PORECHOP } from '../../../modules/nf-core/porechop/porechop/main'
include { PIGZ_UNCOMPRESS   } from '../../../modules/nf-core/pigz/uncompress/main'
include { METAMDBG_ASM      } from '../../../modules/nf-core/metamdbg/asm/main'

workflow ONT_ASSEMBLY {
    take:
    ch_input

    main:

    ch_versions = Channel.empty()

    if (!params.skip_porechop) {
        PORECHOP_PORECHOP(
            ch_input.map { meta, fastq_1, _fastq_2 -> [meta, fastq_1] }
        )
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())

        ch_trimmed_reads = PORECHOP_PORECHOP.out.reads
    }
    else {
        ch_trimmed_reads = ch_input.map { meta, fastq_1, _fastq_2 -> [meta, fastq_1] }
    }

    if (params.use_flye) {
        mode = "--nano-corr"
        FLYE(
            ch_trimmed_reads,
            mode,
        )
        ch_versions = ch_versions.mix(FLYE.out.versions.first())

        PIGZ_UNCOMPRESS(
            FLYE.out.gfa
        )
        ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS.out.versions.first())

        ch_contigs = FLYE.out.fasta
        ch_gfa = PIGZ_UNCOMPRESS.out.file
    }
    else {
        METAMDBG_ASM(
            ch_trimmed_reads,
            "ont",
        )
        ch_versions = ch_versions.mix(METAMDBG_ASM.out.versions.first())

        ch_contigs = METAMDBG_ASM.out.contigs
        ch_gfa = []
    }

    emit:
    reads    = ch_trimmed_reads
    contigs  = ch_contigs
    gfa      = ch_gfa
    versions = ch_versions // channel: [ versions.yml ]
}
