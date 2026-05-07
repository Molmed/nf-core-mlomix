#!/usr/bin/env python3

"""Filter low-variance features and plot variance distributions."""

import argparse
import os
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


def read_features_matrix(file_path: str) -> pd.DataFrame:
    """Read feature matrix with header row (feature_id + sample columns)."""
    features_df = pd.read_csv(file_path, sep=",", header=0, index_col=0)

    if features_df.shape[1] < 1:
        raise ValueError(
            "Expected at least 2 columns in feature file, got "
            f"{features_df.shape[1] + 1}"
        )

    return features_df


def plot_variance_distribution(
    variances: pd.Series,
    threshold: Optional[float],
    title: str,
    output_prefix: str,
    threshold_label: str = "Threshold",
) -> None:
    """Save histogram of feature variances as PNG and SVG."""
    plt.figure(figsize=(10, 6))
    plt.hist(variances, bins=50, color="skyblue", edgecolor="black", alpha=0.7)
    plt.title(title)
    plt.xlabel("Variance", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    if threshold is not None:
        plt.axvline(
            x=threshold,
            color="red",
            linestyle="--",
            label=f"{threshold_label} ({threshold})",
            linewidth=2,
        )
    if threshold is not None:
        plt.legend(fontsize=12)
    plt.grid(axis="y", alpha=0.75)
    plt.tight_layout()
    # plt.savefig(f"{output_prefix}.png", dpi=300)
    plt.savefig(f"{output_prefix}.svg", format="svg", dpi=300)
    plt.close()


def remove_low_variance_features(
    features_df: pd.DataFrame, variance_threshold: float
) -> pd.DataFrame:
    """Keep features (rows/probes) with variance > variance_threshold."""
    print(f"Dataframe shape before filtering: {features_df.shape}")
    variances = features_df.var(axis=1)

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

    filtered_df = features_df.loc[variances > variance_threshold, :]
    print(f"Dataframe shape after filtering: {filtered_df.shape}")

    variances_after = filtered_df.var(axis=1)
    print("Variance statistics after filtering:")
    print(variances_after.describe())

    plot_variance_distribution(
        variances=variances_after,
        threshold=variance_threshold,
        title="Variance distribution after variance filtering",
        output_prefix="variance_distribution_after_filtering",
    )

    return filtered_df


def keep_top_variance_features(
    features_df: pd.DataFrame, keep_top_sites: int
) -> pd.DataFrame:
    """Keep the N rows/probes with highest variance."""
    print(f"Dataframe shape before filtering: {features_df.shape}")
    variances = features_df.var(axis=1)

    print("Variance statistics before filtering:")
    print(variances.describe())

    available_sites = len(variances)
    sites_to_keep = min(keep_top_sites, available_sites)
    if keep_top_sites > available_sites:
        print(
            f"Requested top {keep_top_sites} sites, but only "
            f"{available_sites} are available. Keeping all sites."
        )

    cutoff_variance = variances.nlargest(sites_to_keep).min()
    print(f"Keeping top variance sites: {sites_to_keep}")
    print(f"Variance cutoff at rank {sites_to_keep}: {cutoff_variance}")

    plot_variance_distribution(
        variances=variances,
        threshold=cutoff_variance,
        threshold_label="Top-N cutoff",
        title="Variance distribution before variance filtering",
        output_prefix="variance_distribution_before_filtering",
    )

    top_indices = variances.nlargest(sites_to_keep).index
    filtered_df = features_df.loc[top_indices, :]
    print(f"Dataframe shape after filtering: {filtered_df.shape}")

    variances_after = filtered_df.var(axis=1)
    print("Variance statistics after filtering:")
    print(variances_after.describe())

    plot_variance_distribution(
        variances=variances_after,
        threshold=cutoff_variance,
        threshold_label="Top-N cutoff",
        title="Variance distribution after variance filtering",
        output_prefix="variance_distribution_after_filtering",
    )

    return filtered_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove low-variance features and plot distributions"
    )
    parser.add_argument(
        "--matrix", required=True, help="CSV file of feature values"
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
    parser.add_argument(
        "--keep_top_sites",
        type=int,
        default=None,
        help=(
            "If provided, keep only the N highest-variance features and "
            "ignore --variance_threshold."
        ),
    )
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    if args.variance_threshold < 0:
        raise ValueError("--variance_threshold must be >= 0")
    if args.keep_top_sites is not None and args.keep_top_sites <= 0:
        raise ValueError("--keep_top_sites must be >= 1")

    print("=== Filter By Variance ===")
    print(f"Feature matrix file: {args.matrix}")
    if args.keep_top_sites is not None:
        print(
            f"Keep top variance sites: {args.keep_top_sites} "
            "(overrides variance threshold)"
        )
    else:
        print(f"Variance threshold: {args.variance_threshold}")
    print()

    os.makedirs(args.outdir, exist_ok=True)
    os.chdir(args.outdir)

    features_df = read_features_matrix(args.matrix)
    print(f"Input matrix shape: {features_df.shape}\n")

    if args.keep_top_sites is not None:
        filtered_df = keep_top_variance_features(features_df, args.keep_top_sites)
    else:
        filtered_df = remove_low_variance_features(
            features_df, args.variance_threshold
        )

    output_file = "variance_filtered_betas.csv"
    filtered_df.index.name = None
    filtered_df.to_csv(output_file, header=True)
    print(f"Saved variance-filtered beta values to {output_file}")
    print("Saved variance distribution plots (before/after) as PNG and SVG")

    print("\n=== Filtering complete ===")


if __name__ == "__main__":
    main()
