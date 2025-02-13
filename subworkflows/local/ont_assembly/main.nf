// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

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
    log      = FLYE.out.log
    versions = ch_versions // channel: [ versions.yml ]
}
