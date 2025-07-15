process CALCULATE_PER_TAXON_DIVERSITY {
    tag { meta.id }
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python:3.10.17--f84d2f5369fbe4c7'
        : 'community.wave.seqera.io/library/python:3.10.17--7315eeb528bd3e85'}"

    input:
    tuple val(meta), path(lineage_tsv), path(diversity_metrics_csv)

    output:
    tuple val(meta), path("*.taxon_diversity.tsv"), emit: taxon_diversity_tsv

    script:
    template("per-taxon-diversity.py")

    stub:
    """
    touch ${meta.id}.taxon_diversity.tsv
    """
}
