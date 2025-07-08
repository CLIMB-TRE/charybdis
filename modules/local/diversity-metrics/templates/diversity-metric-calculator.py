#! /usr/bin/env python

import sys
import csv
import math


def calculate_pairwise_difference(
    row, min_cov=0, min_allele_depth=0, min_allele_frequency=0.0
):
    """
    Calculate the nucleotide diversity (π) for a given chromosome position.
    """

    to_check = ("a", "c", "g", "t")

    # Calculate the total number of alleles
    n = sum(
        [
            int(row[field])
            for field in to_check
            if int(row[field]) >= min_allele_depth
            and (int(row[field]) / int(row["cov"])) >= min_allele_frequency
        ]
    )

    if n < min_cov:
        return None

    # Calculate the denominator for the nucleotide diversity formula
    denominator = n * (n - 1)
    # If the denominator is 0, return None
    if denominator == 0 or row["cov"] == "0":
        return None

    # Calculate the allele frequency
    allele_freq = sum(
        [
            (int(row[field]) * (int(row[field]) - 1))
            for field in to_check
            if int(row[field]) >= min_allele_depth
            and (int(row[field]) / int(row["cov"])) >= min_allele_frequency
        ]
    )

    # Calculate the nucleotide diversity (π)
    dl = (denominator - allele_freq) / denominator
    return dl


def calculate_shannon_entropy(
    row, min_cov=0, min_allele_depth=0, min_allele_frequency=0.0
):
    """
    Calculate the Shannon entropy for a given chromosome position.
    """

    to_check = ("a", "c", "g", "t", "ds")

    # Calculate the total number of alleles
    n = sum([int(row[field]) for field in to_check])

    if n < min_cov:
        return None

    # Calculate the denominator for the Shannon entropy formula
    # Hl = -E [pi * log(pi)]
    entropy = -sum(
        (int(row[field]) / n) * math.log(int(row[field]) / n)
        for field in to_check
        if int(row[field]) > 0
        and int(row[field]) >= min_allele_depth
        and (int(row[field]) / n) > min_allele_frequency
    )

    return entropy


def run(args):

    pairwise_difference_arrays = {}
    entropy_arrays = {}

    to_check = ("a", "c", "g", "t", "ds")

    # Load the VCF file
    with open(args.tsv, "r") as tsv_file:
        # Headers -> chrom   pos     ins     cov     a       c       g       t       ds      n

        reader = csv.DictReader(tsv_file, delimiter="\\t")

        for row in reader:

            # TODO support insertions
            if int(row["ins"]) > 0:
                continue

            # Skip if the coverage is 0
            if row["cov"] == "0":
                continue

            pairwise_difference = calculate_pairwise_difference(
                row, args.min_cov, args.min_allele_depth, args.min_allele_frequency
            )
            entropy = calculate_shannon_entropy(
                row, args.min_cov, args.min_allele_depth, args.min_allele_frequency
            )

            pairwise_difference_arrays.setdefault(row["chrom"], [])
            entropy_arrays.setdefault(row["chrom"], [])

            if pairwise_difference is not None:
                pairwise_difference_arrays[row["chrom"]].append(pairwise_difference)

            if entropy is not None:
                entropy_arrays[row["chrom"]].append(entropy)
    with open(f"{args.sample}.diversity_metrics.csv", "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["sample", "chrom", "nt_diversity", "shannon_entropy"],
            lineterminator="\\n",
        )
        writer.writeheader()

        for chrom in pairwise_difference_arrays:

            if not pairwise_difference_arrays[chrom]:
                continue

            pi_diversity = sum(
                x / len(pairwise_difference_arrays[chrom])
                for x in pairwise_difference_arrays[chrom]
            )

            shannon_entropy = sum(
                x / len(entropy_arrays[chrom]) for x in entropy_arrays[chrom]
            )

            out_dict = {
                "sample": args.sample,
                "chrom": chrom,
                "nt_diversity": pi_diversity,
                "shannon_entropy": shannon_entropy,
            }
            writer.writerow(out_dict)


if __name__ == "__main__":
    # import argparse
    from types import SimpleNamespace

    args = SimpleNamespace()
    args.sample = "${meta.id}"
    args.tsv = "${pileup_tsv}"
    args.min_cov = int("${params.diversity_metric_min_cov}")
    args.min_allele_depth = int("${params.diversity_metric_min_allele_depth}")
    args.min_allele_frequency = float("${params.diversity_metric_min_allele_frequency}")

    # parser = argparse.ArgumentParser(
    #     description="Calculate nucleotide diversity from a maptide TSV output file."
    # )
    # parser.add_argument(
    #     "sample", help="Sample name to calculate nucleotide diversity for."
    # )
    # parser.add_argument("tsv", help="Path to the input TSV file (maptide output).")
    # parser.add_argument(
    #     "--min-cov",
    #     type=int,
    #     default=0,
    #     help="Minimum coverage to consider a position.",
    # )
    # parser.add_argument(
    #     "--min-allele-depth",
    #     type=int,
    #     default=0,
    #     help="Minimum allele depth for consideration.",
    # )
    # parser.add_argument(
    #     "--min-allele-frequency",
    #     type=float,
    #     default=0.0,
    #     help="Minimum allele frequency to consider an alternative allele.",
    # )
    # args = parser.parse_args()

    run(args)
