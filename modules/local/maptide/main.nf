process MAPTIDE_PILEUP {
    tag { meta.id }
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pip_maptide:0dd6a062b156b998'
        : 'community.wave.seqera.io/library/pip_maptide:70fcb5df405b3948'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pileup.tsv"), emit: pileup_tsv
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    maptide \\
        ${args} \\
        ${bam} \\
        > ${prefix}.pileup.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        maptide: \$(maptide --version)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileup.tsv
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        maptide: \$(maptide --version)
    END_VERSIONS
    """
}
