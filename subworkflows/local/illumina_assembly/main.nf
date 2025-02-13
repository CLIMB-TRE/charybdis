include { MEGAHIT       } from '../../../modules/nf-core/megahit/main'
include { MEGAHIT_FASTG } from '../../../modules/local/megahit/fastg/main'

workflow ILLUMINA_ASSEMBLY {
    take:
    ch_input

    main:

    ch_versions = Channel.empty()

    ch_fastq = ch_input.branch { _id, platform, _fastq_1, _fastq_2 ->
        paired: platform == "illumina"
        single: platform == "illumina.se"
    }

    ch_fastq.paired
        .map { meta, platform, fastq_1, fastq_2 -> [[id: meta.id, platform: platform, single_end: false], fastq_1, fastq_2] }
        .set { ch_paired }

    ch_fastq.single
        .map { meta, platform, fastq_1, _fastq_2 -> [[id: meta.id, platform: platform, single_end: true], fastq_1, [:]] }
        .set { ch_single }

    ch_to_assemble = ch_paired.mix(ch_single)

    MEGAHIT(ch_to_assemble)
    ch_versions = ch_versions.mix(MEGAHIT.versions.first())

    MEGAHIT_FASTG(MEGAHIT.out.kfinal_contigs)
    ch_versions = ch_versions.mix(MEGAHIT_FASTG.versions.first())

    emit:
    contigs  = MEGAHIT.out.contigs
    fastg    = MEGAHIT_FASTG.out.fastg
    versions = ch_versions // channel: [ path(versions.yml) ]
}
