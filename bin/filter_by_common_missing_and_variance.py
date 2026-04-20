#!/usr/bin/env python3

"""Apply common-probe, missingness, and variance filtering in one pass."""

import argparse
import os

import matplotlib.pyplot as plt  # type: ignore[reportMissingImports]
import pandas as pd


def deduplicate_probe_columns(beta_df: pd.DataFrame) -> pd.DataFrame:
    """Average duplicate columns (e.g. EPIC v2 replicate probes)."""
    n_before = beta_df.shape[1]
    beta_df = beta_df.groupby(beta_df.columns, axis=1).mean()
    n_after = beta_df.shape[1]
    n_dupes = n_before - n_after
    if n_dupes > 0:
        print(
            f"Deduplicated {n_dupes} duplicate probe columns "
            f"({n_before} -> {n_after})"
        )
    else:
        print("No duplicate probe columns found")
    return beta_df


def filter_to_common_probes(
    beta_df: pd.DataFrame,
    probe_list: list[str],
) -> pd.DataFrame:
    """Filter rows to probes present in the provided reference list."""
    probe_intersection = beta_df.index.intersection(probe_list)
    n_before = beta_df.shape[0]
    beta_df = beta_df.loc[probe_intersection]
    n_after = beta_df.shape[0]
    print(
        f"Filtered probes: {n_before} -> {n_after} "
        f"({n_before - n_after} removed)"
    )
    return beta_df


def remove_high_missing_features(
    beta_df: pd.DataFrame,
    missing_threshold: float,
) -> pd.DataFrame:
    """Keep columns with missing fraction <= missing_threshold."""
    print(f"Dataframe shape before missingness filtering: {beta_df.shape}")
    percentage_missing_vals = beta_df.isnull().mean()
    print("Missingness statistics (percentage of missing values per feature):")
    print(percentage_missing_vals.describe())
    filtered_df = beta_df.loc[:, percentage_missing_vals <= missing_threshold]
    print(f"Dataframe shape after missingness filtering: {filtered_df.shape}")
    return filtered_df


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
    beta_df: pd.DataFrame,
    variance_threshold: float,
) -> pd.DataFrame:
    """Keep columns with variance > variance_threshold."""
    print(f"Dataframe shape before variance filtering: {beta_df.shape}")
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
    print(f"Dataframe shape after variance filtering: {filtered_df.shape}")

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
        description=(
            "Run common-probe, missingness, and variance filters "
            "in sequence"
        )
    )
    parser.add_argument(
        "--betas",
        required=True,
        help="TSV file of beta values",
    )
    parser.add_argument(
        "--probes",
        required=True,
        help="Text file with one probe ID per line",
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

    if not 0 <= args.missing_threshold <= 1:
        raise ValueError("--missing_threshold must be between 0 and 1")
    if args.variance_threshold < 0:
        raise ValueError("--variance_threshold must be >= 0")

    print("=== Combined DNAm Filtering ===")
    print(f"Beta values file: {args.betas}")
    print(f"Probe list file: {args.probes}")
    print(f"Missingness threshold: {args.missing_threshold}")
    print(f"Variance threshold: {args.variance_threshold}\n")

    os.makedirs(args.outdir, exist_ok=True)
    os.chdir(args.outdir)

    beta_df = pd.read_csv(args.betas, sep="\t", index_col=0)
    print(f"Input matrix shape: {beta_df.shape}\n")

    print("Step 1: Deduplicate probes...")
    beta_df = deduplicate_probe_columns(beta_df)
    print(f"Matrix after deduplication: {beta_df.shape}\n")

    print("Step 2: Filter to common probes...")
    probes_to_keep = pd.read_csv(
        args.probes,
        sep="\t",
        header=None,
        names=["probe_id"],
    )
    probe_list = probes_to_keep["probe_id"].tolist()
    print(f"Reference probe list: {len(probe_list)} probes")
    beta_df = filter_to_common_probes(beta_df, probe_list)
    print(f"Matrix after common probe filtering: {beta_df.shape}\n")

    print("Step 3: Filter by missingness...")
    beta_df = remove_high_missing_features(beta_df, args.missing_threshold)
    print(f"Matrix after missingness filtering: {beta_df.shape}\n")

    print("Step 4: Filter by variance and plot distributions...")
    beta_df = remove_low_variance_features(beta_df, args.variance_threshold)
    print(f"Matrix after variance filtering: {beta_df.shape}\n")

    output_file = "variance_filtered_betas.tsv"
    beta_df.to_csv(
        output_file,
        sep="\t",
        index_label=beta_df.index.name or "probe_id",
    )
    print(f"Saved final filtered beta values to {output_file}")
    print("Saved variance distribution plots (before/after) as PNG and SVG")

    print("\n=== Filtering complete ===")


if __name__ == "__main__":
    main()
