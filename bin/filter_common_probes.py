#!/usr/bin/env python3

"""
Filter common probes across methylation array platforms.

1. Deduplicates probes by averaging duplicate columns (e.g. EPIC v2 replicate probes)
2. Filters to only probes present in a reference probe list (e.g. probes common to 450K, EPIC v1, EPIC v2)
"""

import argparse
import os

import numpy as np
import pandas as pd


def deduplicate_probes(beta_df: pd.DataFrame) -> pd.DataFrame:
    """Average duplicate probe columns (e.g. EPIC v2 replicate probes)."""
    n_before = beta_df.shape[1]
    beta_df = beta_df.groupby(beta_df.columns, axis=1).mean()
    n_after = beta_df.shape[1]
    n_dupes = n_before - n_after
    if n_dupes > 0:
        print(f"Deduplicated {n_dupes} duplicate probe columns ({n_before} -> {n_after})")
    else:
        print("No duplicate probe columns found")
    return beta_df


def filter_probes(beta_df: pd.DataFrame, probe_list: list) -> pd.DataFrame:
    """Filter beta matrix to only probes in the reference list."""
    probe_intersection = beta_df.index.intersection(probe_list)
    n_before = beta_df.shape[0]
    beta_df = beta_df.loc[probe_intersection]
    n_after = beta_df.shape[0]
    print(f"Filtered probes: {n_before} -> {n_after} ({n_before - n_after} removed)")
    return beta_df


def main():
    parser = argparse.ArgumentParser(description="Filter common probes across platforms")
    parser.add_argument("--betas", required=True, help="TSV file of beta values (probe_id as index)")
    parser.add_argument("--probes", required=True, help="Text file with one probe ID per line")
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    print("=== Filter Common Probes ===")
    print(f"Beta values file: {args.betas}")
    print(f"Probe list file: {args.probes}\n")

    # Read beta values
    beta_df = pd.read_csv(args.betas, sep="\t", index_col=0)
    print(f"Input beta matrix: {beta_df.shape[0]} probes x {beta_df.shape[1]} samples\n")

    # Step 1: Deduplicate probes
    print("Step 1: Deduplicating probes...")
    beta_df = deduplicate_probes(beta_df)
    print(f"Matrix after dedup: {beta_df.shape[0]} probes x {beta_df.shape[1]} samples\n")

    # Step 2: Filter to common probes
    print("Step 2: Filtering to common probes...")
    probes_to_keep = pd.read_csv(args.probes, sep="\t", header=None, names=["probe_id"])
    probe_list = probes_to_keep["probe_id"].tolist()
    print(f"Reference probe list: {len(probe_list)} probes")
    beta_df = filter_probes(beta_df, probe_list)
    print(f"Matrix after filtering: {beta_df.shape[0]} probes x {beta_df.shape[1]} samples\n")

    # Save output
    output_file = os.path.join(args.outdir, "filtered_betas.tsv")
    beta_df.to_csv(output_file, sep="\t", index_label="probe_id")
    print(f"Saved filtered beta values to {output_file}")

    print("\n=== Filtering complete ===")


if __name__ == "__main__":
    main()
