process DIVERSITY_METRIC_PLOT {
    tag { meta.id }
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/numpy_pandas_plotly:ed8e09205d7f02b8'
        : 'community.wave.seqera.io/library/numpy_pandas_plotly:9dad17f2461b50ef'}"

    errorStrategy { task.exitStatus == 15 ? 'ignore' : 'terminate' }

    input:
    tuple val(meta), path(diversity_metric_tsv)

    output:
    tuple val(meta), path("*.species_diversity_metric_plot.html"), emit: diversity_metric_plot_species_html
    tuple val(meta), path("*.genus_diversity_metric_plot.html"), emit: diversity_metric_plot_genus_html

    script:
    template("diversity-metric-plot.py")

    stub:
    """
    touch ${meta.id}.species_diversity_metric_plot.html
    touch ${meta.id}.genus_diversity_metric_plot.html
    """
}
