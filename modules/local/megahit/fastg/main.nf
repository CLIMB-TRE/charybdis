process MEGAHIT_FASTG {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f2cb827988dca7067ff8096c37cb20bc841c878013da52ad47a50865d54efe83/data'
        : 'community.wave.seqera.io/library/megahit_pigz:87a590163e594224'}"

    input:
    tuple val(meta), path(kfinal_contigs)

    output:
    tuple val(meta), path("*.fastg"), emit: fastg
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    megahit_toolkit contig2fastg \\
        ${params.bandage_kmer_size} \\
        k${params.bandage_kmer_size}.contigs.fa \\
        ${args} \\
        > ${prefix}.fastg

    if [ ! -s ${prefix}.fastg ]; then
        echo "ERROR: Failed to generate FASTG file (probably due to insufficient contig overlaps)"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit_toolkit: \$(megahit_toolkit dumpversion)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit_toolkit: \$(megahit_toolkit dumpversion)
    END_VERSIONS
    """
}
