process DIVERSITY_METRICS {
    tag { meta.id }
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python:3.10.17--f84d2f5369fbe4c7'
        : 'community.wave.seqera.io/library/python:3.10.17--7315eeb528bd3e85'}"

    input:
    tuple val(meta), path(pileup_tsv)

    output:
    tuple val(meta), path("${meta.id}.diversity_metrics.csv"), emit: csv

    script:
    template("diversity-metric-calculator.py")

    stub:
    """
    touch ${meta.id}.diversity_metrics.csv
    """
}
