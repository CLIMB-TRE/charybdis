include { FLYE                                       } from '../../../modules/nf-core/flye/main'
include { PORECHOP_PORECHOP                          } from '../../../modules/nf-core/porechop/porechop/main'
include { PIGZ_UNCOMPRESS as PIGZ_UNCOMPRESS_FLYE    } from '../../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as PIGZ_UNCOMPRESS_MINIASM } from '../../../modules/nf-core/pigz/uncompress/main'
include { METAMDBG_ASM                               } from '../../../modules/nf-core/metamdbg/asm/main'
include { MINIASM                                    } from '../../../modules/nf-core/miniasm/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_MINIASM   } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_RACON     } from '../../../modules/nf-core/minimap2/align/main'
include { RACON                                      } from '../../../modules/nf-core/racon/main'

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

    if (params.ont_assembler == "flye") {
        mode = "--nano-corr"
        FLYE(
            ch_trimmed_reads,
            mode,
        )
        ch_versions = ch_versions.mix(FLYE.out.versions.first())

        PIGZ_UNCOMPRESS_FLYE(
            FLYE.out.gfa
        )
        ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS_FLYE.out.versions.first())

        ch_contigs = FLYE.out.fasta
        ch_gfa = PIGZ_UNCOMPRESS_FLYE.out.file
    }
    else if (params.ont_assembler == "metamdbg") {
        METAMDBG_ASM(
            ch_trimmed_reads,
            "ont",
        )
        ch_versions = ch_versions.mix(METAMDBG_ASM.out.versions.first())

        ch_contigs = METAMDBG_ASM.out.contigs
        ch_gfa = []
    }
    else if (params.ont_assembler == "miniasm") {
        MINIMAP2_ALIGN_MINIASM(
            ch_trimmed_reads,
            [[:], []],
            false,
            false,
            false,
            false,
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_MINIASM.out.versions.first())

        ch_miniasm = ch_trimmed_reads.join(MINIMAP2_ALIGN_MINIASM.out.paf, by: 0)

        MINIASM(
            ch_miniasm
        )
        ch_versions = ch_versions.mix(MINIASM.out.versions.first())

        PIGZ_UNCOMPRESS_MINIASM(
            MINIASM.out.gfa
        )
        ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS_MINIASM.out.versions.first())

        MINIMAP2_ALIGN_RACON(
            ch_trimmed_reads,
            MINIASM.out.assembly,
            false,
            false,
            false,
            false,
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_RACON.out.versions.first())

        ch_racon = ch_trimmed_reads.join(MINIASM.out.assembly, by: 0).join(MINIMAP2_ALIGN_RACON.out.paf, by: 0)

        RACON(ch_racon)

        ch_contigs = RACON.out.improved_assembly
        ch_gfa = PIGZ_UNCOMPRESS_MINIASM.out.file
    }
    else {
        error("Unrecognised assembler: ${params.ont_assembler}")
    }

    emit:
    reads    = ch_trimmed_reads
    contigs  = ch_contigs
    gfa      = ch_gfa
    versions = ch_versions // channel: [ versions.yml ]
}
