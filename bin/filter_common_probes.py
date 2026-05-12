#!/usr/bin/env python3

"""
Filter common probes across methylation array platforms.

1. Filters to probes present in a reference probe list
   (e.g. probes common to 450K, EPIC v1, EPIC v2).
"""

import argparse
import os
import re

import pandas as pd


def normalize_probe_id(probe_id: str) -> str:
    """Strip any underscore suffix and trailing numeric version suffix."""
    if isinstance(probe_id, str):
        probe_id = probe_id.split("_", 1)[0]
        return re.sub(r"\.\d+$", "", probe_id)
    return probe_id


def read_beta_matrix(
    file_path: str,
    sample_name: str | None = None,
) -> pd.DataFrame:
    """Read headerless beta matrix (probe_id + beta columns)."""
    beta_df = pd.read_csv(file_path, sep="\t", header=None)

    if beta_df.shape[1] < 2:
        raise ValueError(
            "Expected at least 2 columns in beta file, got "
            f"{beta_df.shape[1]}"
        )

    beta_df.columns = ["probe_id"] + [
        f"value_{i}" for i in range(1, beta_df.shape[1])
    ]
    if sample_name and beta_df.shape[1] == 2:
        beta_df = beta_df.rename(columns={"value_1": sample_name})

    beta_df = beta_df.set_index("probe_id")
    beta_df.index = beta_df.index.map(normalize_probe_id)
    return beta_df


def filter_probes(beta_df: pd.DataFrame, probe_list: list) -> pd.DataFrame:
    """Filter beta matrix to only probes in the reference list."""
    probe_intersection = beta_df.index.intersection(probe_list)
    n_before = beta_df.shape[0]
    beta_df = beta_df.loc[probe_intersection]
    n_after = beta_df.shape[0]
    print(
        f"Filtered probes: {n_before} -> {n_after} "
        f"({n_before - n_after} removed)"
    )
    return beta_df


def main():
    parser = argparse.ArgumentParser(
        description="Filter common probes across platforms"
    )
    parser.add_argument(
        "--betas",
        required=True,
        help="TSV file of beta values (probe_id as index)",
    )
    parser.add_argument(
        "--probes",
        required=True,
        help="Text file with one probe ID per line",
    )
    parser.add_argument(
        "--sample-name",
        default=None,
        help="Optional sample name for single-sample input files",
    )
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    print("=== Filter Common Probes ===")
    print(f"Beta values file: {args.betas}")
    print(f"Probe list file: {args.probes}\n")

    # Read beta values
    beta_df = read_beta_matrix(args.betas, args.sample_name)
    print(
        f"Input beta matrix: {beta_df.shape[0]} probes x "
        f"{beta_df.shape[1]} samples\n"
    )

    # Filter to common probes
    print("Filtering to common probes...")
    probes_to_keep = pd.read_csv(
        args.probes,
        sep="\t",
        header=None,
        names=["probe_id"],
    )
    probe_list = probes_to_keep["probe_id"].map(normalize_probe_id).tolist()
    print(f"Reference probe list: {len(probe_list)} probes")
    beta_df = filter_probes(beta_df, probe_list)
    print(
        f"Matrix after filtering: {beta_df.shape[0]} probes x "
        f"{beta_df.shape[1]} samples\n"
    )

    # Save output
    output_file = os.path.join(args.outdir, "filtered_betas.tsv")
    beta_df.index.name = None
    beta_df.to_csv(output_file, sep="\t", header=False)
    print(f"Saved filtered beta values to {output_file}")

    print("\n=== Filtering complete ===")


if __name__ == "__main__":
    main()
