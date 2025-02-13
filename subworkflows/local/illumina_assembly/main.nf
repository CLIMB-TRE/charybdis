include { MEGAHIT       } from '../../../modules/nf-core/megahit/main'
include { MEGAHIT_FASTG } from '../../../modules/local/megahit/fastg/main'

workflow ILLUMINA_ASSEMBLY {
    take:
    ch_input

    main:

    ch_versions = Channel.empty()

    MEGAHIT(ch_input)
    ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())

    MEGAHIT_FASTG(MEGAHIT.out.k_contigs)
    ch_versions = ch_versions.mix(MEGAHIT_FASTG.out.versions.first())

    emit:
    contigs  = MEGAHIT.out.contigs
    fastg    = MEGAHIT_FASTG.out.fastg
    versions = ch_versions // channel: [ path(versions.yml) ]
}
