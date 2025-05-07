process KRAKEN2_CLIENT {
    tag "${meta.id}"

    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "quay.io/climb-tre/kraken2-server:sha-13ee499"

    input:
    tuple val(meta), path(fastx)
    val k2_remote

    output:
    tuple val(meta), path("${meta.id}.kraken2.classifiedreads.txt"), emit: assignments
    tuple val(meta), path("${meta.id}.kraken2.report.txt"), emit: report
    path ("versions.yml"), emit: versions

    script:

    k2_remote_splits = k2_remote.split(":")

    """
    kraken2_client \
        --host-ip ${k2_remote_splits[0]} --port ${k2_remote_splits[1]} \
        --report "${meta.id}.kraken2.report.txt" \
        --sequence ${fastx} > "${meta.id}.kraken2.classifiedreads.txt"

    # I know this is awful but kraken2_client doesn't output its version for some mysterious reason
    echo "kraken2_client: 0.1.7" > versions.yml
    """
}
