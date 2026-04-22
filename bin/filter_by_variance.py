#!/usr/bin/env python3

"""Filter low-variance features and plot variance distributions."""

import argparse
import os

import matplotlib.pyplot as plt
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


def plot_variance_distribution(
    variances: pd.Series,
    threshold: float,
    title: str,
    output_prefix: str,
) -> None:
    """Save histogram of feature variances as PNG and SVG."""
    plt.figure(figsize=(10, 6))
    plt.hist(variances, bins=50, color="skyblue", edgecolor="black", alpha=0.7)
    plt.title(title)
    plt.xlabel("Variance", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.axvline(
        x=threshold,
        color="red",
        linestyle="--",
        label=f"Threshold ({threshold})",
        linewidth=2,
    )
    plt.legend(fontsize=12)
    plt.grid(axis="y", alpha=0.75)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=300)
    plt.savefig(f"{output_prefix}.svg", format="svg", dpi=300)
    plt.close()


def remove_low_variance_features(
    beta_df: pd.DataFrame, variance_threshold: float
) -> pd.DataFrame:
    """Keep features (columns) with variance > variance_threshold."""
    print(f"Dataframe shape before filtering: {beta_df.shape}")
    variances = beta_df.var(axis=0)

    print("Variance statistics before filtering:")
    print(variances.describe())

    num_features_above = (variances > variance_threshold).sum()
    frequency_density = (variances > variance_threshold).mean()
    print(
        f"Number of features above {variance_threshold}: "
        f"{num_features_above}"
    )
    print(
        f"Frequency density for variance > {variance_threshold}: "
        f"{frequency_density:.2%}"
    )

    plot_variance_distribution(
        variances=variances,
        threshold=variance_threshold,
        title="Variance distribution before variance filtering",
        output_prefix="variance_distribution_before_filtering",
    )

    filtered_df = beta_df.loc[:, variances > variance_threshold]
    print(f"Dataframe shape after filtering: {filtered_df.shape}")

    variances_after = filtered_df.var(axis=0)
    print("Variance statistics after filtering:")
    print(variances_after.describe())

    plot_variance_distribution(
        variances=variances_after,
        threshold=variance_threshold,
        title="Variance distribution after variance filtering",
        output_prefix="variance_distribution_after_filtering",
    )

    return filtered_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove low-variance features and plot distributions"
    )
    parser.add_argument(
        "--betas", required=True, help="TSV file of beta values"
    )
    parser.add_argument(
        "--variance_threshold",
        type=float,
        default=0.01,
        help=(
            "Minimum variance required for features to be retained "
            "(default: 0.01)"
        ),
    )
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    if args.variance_threshold < 0:
        raise ValueError("--variance_threshold must be >= 0")

    print("=== Filter By Variance ===")
    print(f"Beta values file: {args.betas}")
    print(f"Variance threshold: {args.variance_threshold}\n")

    os.makedirs(args.outdir, exist_ok=True)
    os.chdir(args.outdir)

    beta_df = read_beta_matrix(args.betas)
    print(f"Input matrix shape: {beta_df.shape}\n")

    filtered_df = remove_low_variance_features(
        beta_df, args.variance_threshold
    )

    output_file = "variance_filtered_betas.tsv"
    filtered_df.index.name = None
    filtered_df.to_csv(output_file, sep="\t", header=True)
    print(f"Saved variance-filtered beta values to {output_file}")
    print("Saved variance distribution plots (before/after) as PNG and SVG")

    print("\n=== Filtering complete ===")


if __name__ == "__main__":
    main()
