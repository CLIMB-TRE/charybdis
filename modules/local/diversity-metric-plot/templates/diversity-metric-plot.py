#! /usr/bin/env python

import csv
import plotly.express as px
import pandas as pd
import sys


def parse_metrics_tsv(metrics_tsv: str) -> dict:
    """
    Parse the metrics TSV file and return a dictionary mapping taxon IDs to their nucleotide diversity scores.

    Args:
        metrics_tsv (str): Path to the metrics TSV file.

    Returns:
        list: A list of dictionaries, each containing taxon ID, diversity score, and associated info.
    """
    taxon_scores = {}

    with open(metrics_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        taxon_scores = [row for row in reader]

    return taxon_scores


def diversity_metric_plot(
    taxon_scores: pd.DataFrame, output_file: str, rank: str = "species"
):
    """
    Generate a diversity metric plot from the taxon scores.

    Args:
        taxon_scores (pd.DataFrame): DataFrame containing taxon scores.
        output_file (str): Path to the output file for the plot.
    """

    filtered_taxon_scores = taxon_scores[taxon_scores["rank"] == rank]

    fig = px.scatter(
        filtered_taxon_scores,
        x="mean_depth_zscore",
        y="nt_diversity_zscore",
        color="scientific_name",
        hover_name="tax_id",
        title="Abundance / Diversity Metric Plot",
        labels={
            "nt_diversity_zscore": "Nucleotide Diversity Z-Score",
            "mean_depth_zscore": "Mean Depth Z-Score",
        },
        color_discrete_sequence=px.colors.qualitative.Plotly,
    )

    fig.update_layout(
        xaxis_title="Mean Depth Z-Score",
        yaxis_title="Nucleotide Diversity Z-Score",
        legend_title="Scientific Name",
        template="plotly_white",
        yaxis=dict(type="linear", categoryorder="total ascending"),
        xaxis=dict(type="linear", categoryorder="total ascending"),
    )

    # fig = px.scatter(
    #     taxon_scores,
    #     x="rank",
    #     y="relative_nt_diversity",
    #     color="scientific_name",
    #     hover_name="tax_id",
    #     title="Diversity Metric Plot",
    #     labels={
    #         "relative_nt_diversity": "Relative Nucleotide Diversity",
    #         "rank": "Taxonomic Rank",
    #     },
    #     color_discrete_sequence=px.colors.qualitative.Plotly,
    # )

    # fig.update_layout(
    #     yaxis=dict(type="linear", categoryorder="total ascending"),
    #     xaxis=dict(
    #         type="category",
    #         categoryarray=[
    #             "no rank",
    #             "superkingdom",
    #             "kingdom",
    #             "clade",
    #             "subkingdom",
    #             "phylum",
    #             "subphylum",
    #             "superclass",
    #             "class",
    #             "subclass",
    #             "superorder",
    #             "infraorder",
    #             "order",
    #             "suborder",
    #             "parvorder",
    #             "infraorder",
    #             "superfamily",
    #             "family",
    #             "subfamily",
    #             "genus",
    #             "species group",
    #             "species complex",
    #             "species",
    #             "subspecies",
    #             "strain",
    #         ],
    #     ),
    # )

    # fig.update_layout(
    #     xaxis_title="Taxonomic Rank",
    #     yaxis_title="Relative Nucleotide Diversity",
    #     legend_title="Scientific Name",
    #     template="plotly_white",
    # )

    fig.write_html(output_file)


def main():
    taxon_scores = parse_metrics_tsv("${diversity_metric_tsv}")
    taxon_score_df = pd.DataFrame(taxon_scores)

    if len(taxon_score_df) == 0:
        print("No taxon scores found in the provided TSV file, cannot generate plot.")
        sys.exit(15)

    diversity_metric_plot(
        taxon_score_df, "${meta.id}.species_diversity_metric_plot.html", rank="species"
    )

    diversity_metric_plot(
        taxon_score_df, "${meta.id}.genus_diversity_metric_plot.html", rank="genus"
    )


if __name__ == "__main__":
    main()
