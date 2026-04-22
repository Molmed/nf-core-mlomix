#!/usr/bin/env python3

"""Filter features by missingness in a beta-value matrix."""

import argparse
import os

import pandas as pd


def read_beta_matrix(file_path: str) -> pd.DataFrame:
    """Read beta matrix with header row (probe_id + sample columns)."""
    beta_df = pd.read_csv(file_path, sep="\t", header=0, index_col=0)

    if beta_df.shape[1] < 1:
        raise ValueError(
            "Expected at least 2 columns in beta file, got "
            f"{beta_df.shape[1] + 1}"
        )

    return beta_df


def remove_high_missing_features(
    beta_df: pd.DataFrame, missing_threshold: float
) -> pd.DataFrame:
    """Keep features (columns) with missing fraction <= missing_threshold."""
    print(f"Dataframe shape before filtering: {beta_df.shape}")
    percentage_missing_vals = beta_df.isnull().mean()
    # Log percentage_missing_vals.describe() for debugging
    print("Missingness statistics (percentage of missing values per feature):")
    print(percentage_missing_vals.describe())
    filtered_df = beta_df.loc[:, percentage_missing_vals <= missing_threshold]
    print(f"Dataframe shape after filtering: {filtered_df.shape}")
    return filtered_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove features with too many missing values"
    )
    parser.add_argument(
        "--betas", required=True, help="TSV file of beta values"
    )
    parser.add_argument(
        "--missing_threshold",
        type=float,
        default=0.1,
        help=(
            "Maximum allowed fraction of missing values per feature "
            "(default: 0.1)"
        ),
    )
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    if not 0 <= args.missing_threshold <= 1:
        raise ValueError("--missing_threshold must be between 0 and 1")

    print("=== Filter By Missingness ===")
    print(f"Beta values file: {args.betas}")
    print(f"Missingness threshold: {args.missing_threshold}\n")

    beta_df = read_beta_matrix(args.betas)
    print(f"Input matrix shape: {beta_df.shape}\n")

    filtered_df = remove_high_missing_features(beta_df, args.missing_threshold)

    output_file = os.path.join(args.outdir, "missing_filtered_betas.tsv")
    filtered_df.index.name = None
    filtered_df.to_csv(output_file, sep="\t", header=True)
    print(f"Saved missingness-filtered beta values to {output_file}")

    print("\n=== Filtering complete ===")


if __name__ == "__main__":
    main()
