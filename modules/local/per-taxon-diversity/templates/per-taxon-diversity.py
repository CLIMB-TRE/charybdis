#! /usr/bin/env python
import csv
import sys
import statistics


def parse_lineage_tsv(lineage_tsv: str) -> dict:
    """
    Parse the lineage TSV file and return a dictionary mapping contig IDs to lineages / a set of taxIDs.

    Args:
        lineage_tsv (str): Path to the lineage TSV file.

    Returns:
        tuple: A tuple containing:
            - taxon_dict (dict): A dictionary mapping taxon IDs to their scientific names,
                                 ranks, and nucleotide diversity scores (fully populated downstream)
            - contig_lookup (dict): A dictionary mapping contig IDs to their lineage tax IDs

    Raises:
        ValueError: If the lineage columns are not of equal length.
    """
    taxon_dict = {}
    contig_lookup = {}

    with open(lineage_tsv, "r") as f:
        csv.field_size_limit(sys.maxsize)
        reader = csv.reader(f, delimiter="\\t")
        for row in reader:
            if row[0].startswith("#") or row[0].startswith("U"):
                continue

            contig_id = row[1]
            tax_id = int(row[2])
            if not all(row[5:8]):
                # Skip rows with missing lineage information
                print(
                    f"Skipping contig {contig_id} due to missing lineage information.",
                    file=sys.stderr,
                )
                continue

            if int(row[3]) < int("${params.contig_min_length}"):
                print(
                    f"Skipping contig {contig_id} due to insufficient length.",
                    file=sys.stderr,
                )
                continue

            lineage_scientific_name = row[5].split(";")
            lineage_tax_id = [int(x) for x in row[6].split(";")]
            lineage_ranks = row[7].split(";")

            if (
                len(lineage_scientific_name)
                != len(lineage_tax_id)
                != len(lineage_ranks)
            ):
                raise ValueError("Lineage columns are not of equal length")

            contig_lookup[contig_id] = lineage_tax_id
            for i, tax_id in enumerate(lineage_tax_id):
                taxon_dict.setdefault(
                    tax_id,
                    {
                        "scientific_name": lineage_scientific_name[i],
                        "rank": lineage_ranks[i],
                        "nt_diversity_scores": [],
                        "contigs": [],
                        "mean_depths": [],
                        "contributing_contigs": 0,
                    },
                )

    return taxon_dict, contig_lookup


def calculate_per_taxon_diversity(
    taxon_dict: dict,
    contig_lookup: dict,
    diversity_metrics_csv: str,
) -> dict:
    """
    Calculate per-taxon nucleotide diversity from the diversity metrics CSV file.

    Args:
        taxon_dict (dict): A dictionary mapping taxon IDs to their scientific names, ranks
                        and nucleotide diversity scores.
        contig_lookup (dict): A dictionary mapping contig IDs to their lineage tax IDs.
        diversity_metrics_csv (str): Path to the diversity metrics CSV file.

    Returns:
        dict: A dictionary mapping taxon IDs to their average nucleotide diversity scores,
                contigs, and the number of contributing contigs.
    """

    with open(diversity_metrics_csv, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            contig_id = row["chrom"]
            nt_diversity = float(row["nt_diversity"])
            tax_ids = contig_lookup.get(contig_id, [])

            for tax_id in tax_ids:
                if tax_id not in taxon_dict:
                    continue

                taxon_dict[tax_id]["nt_diversity_scores"].append(nt_diversity)
                taxon_dict[tax_id]["contigs"].append(contig_id)
                taxon_dict[tax_id]["mean_depths"].append(float(row["mean_depth"]))
                taxon_dict[tax_id]["contributing_contigs"] += 1

    # Calculate average nucleotide diversity for each taxon
    for tax_id, taxon_info in taxon_dict.items():
        if taxon_info["nt_diversity_scores"]:
            avg_diversity = sum(taxon_info["nt_diversity_scores"]) / len(
                taxon_info["nt_diversity_scores"]
            )
        else:
            avg_diversity = 0.0

        taxon_info["nt_diversity"] = avg_diversity

        if taxon_info["mean_depths"]:
            avg_depth = sum(taxon_info["mean_depths"]) / len(taxon_info["mean_depths"])
        else:
            avg_depth = 0.0

        taxon_info["mean_depth"] = avg_depth

    # Calculate nt diversity Z score per taxon (std deviations from mean)
    mean_diversity = (
        sum(taxon_info["nt_diversity"] for taxon_info in taxon_dict.values())
        / len(taxon_dict)
        if taxon_dict
        else 0.0
    )

    std_dev_diversity = (
        statistics.stdev(
            taxon_info["nt_diversity"] for taxon_info in taxon_dict.values()
        )
        if len(taxon_dict) > 1
        else 0.0
    )

    for taxon_info in taxon_dict.values():
        if std_dev_diversity > 0:
            taxon_info["nt_diversity_zscore"] = (
                taxon_info["nt_diversity"] - mean_diversity
            ) / std_dev_diversity
        else:
            taxon_info["nt_diversity_zscore"] = ""

    mean_depth = (
        sum(taxon_info["mean_depth"] for taxon_info in taxon_dict.values())
        / len(taxon_dict)
        if taxon_dict
        else 0.0
    )

    std_dev_depth = (
        statistics.stdev(taxon_info["mean_depth"] for taxon_info in taxon_dict.values())
        if len(taxon_dict) > 1
        else 0.0
    )

    for taxon_info in taxon_dict.values():
        if std_dev_depth > 0:
            taxon_info["mean_depth_zscore"] = (
                taxon_info["mean_depth"] - mean_depth
            ) / std_dev_depth
        else:
            taxon_info["mean_depth_zscore"] = ""

    return taxon_dict


def write_per_taxon_diversity(
    taxon_dict: dict,
    output_file: str,
):
    """Write the per-taxon diversity results to a TSV file."""
    with open(output_file, "w") as f:
        writer = csv.DictWriter(
            f,
            delimiter="\\t",
            lineterminator="\\n",
            fieldnames=[
                "tax_id",
                "scientific_name",
                "rank",
                "nt_diversity",
                "nt_diversity_zscore",
                "mean_contig_depth",
                "mean_depth_zscore",
                "contributing_contigs",
            ],
        )

        writer.writeheader()

        out_rows = []

        for tax_id, taxon_info in taxon_dict.items():
            print(taxon_info)
            out_rows.append(
                {
                    "tax_id": tax_id,
                    "scientific_name": taxon_info["scientific_name"],
                    "rank": taxon_info["rank"],
                    "nt_diversity": taxon_info["nt_diversity"],
                    "nt_diversity_zscore": taxon_info["nt_diversity_zscore"],
                    "mean_contig_depth": taxon_info["mean_depth"],
                    "mean_depth_zscore": taxon_info["mean_depth_zscore"],
                    "contributing_contigs": taxon_info["contributing_contigs"],
                }
            )

        out_rows.sort(key=lambda x: x["tax_id"])

        for row in out_rows:
            writer.writerow(row)


def main():
    import argparse

    from types import SimpleNamespace

    args = SimpleNamespace()
    args.lineage_tsv = "${lineage_tsv}"
    args.diversity_metrics_csv = "${diversity_metrics_csv}"
    args.output_file = "${meta.id}.taxon_diversity.tsv"

    # parser = argparse.ArgumentParser(
    #     description="Calculate per-taxon nucleotide diversity from a maptide TSV output file."
    # )
    # parser.add_argument(
    #     "--lineage_tsv",
    #     required=True,
    #     help="Path to the lineage TSV file.",
    # )
    # parser.add_argument(
    #     "--diversity_metrics_csv",
    #     required=True,
    #     help="Path to the diversity metrics CSV file.",
    # )
    # parser.add_argument(
    #     "--output_file",
    #     required=True,
    #     help="Output file for per-taxon diversity results.",
    # )

    # args = parser.parse_args()

    taxon_dict, contig_lookup = parse_lineage_tsv(args.lineage_tsv)
    print(
        f"Parsed {len(taxon_dict)} taxa and {len(contig_lookup)} contigs from lineage TSV.",
        file=sys.stderr,
    )
    taxon_dict = calculate_per_taxon_diversity(
        taxon_dict, contig_lookup, args.diversity_metrics_csv
    )
    print(
        f"Calculated per-taxon diversity for {len(taxon_dict)} taxa.",
        file=sys.stderr,
    )
    write_per_taxon_diversity(taxon_dict, args.output_file)
    print(
        f"Wrote per-taxon diversity results to {args.output_file}.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
