#!/usr/bin/env python3

"""
P-value correction for methylation beta values.
Replaces beta values with NaN where detection p-value >= threshold.
"""

import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="P-value correction for beta values"
    )
    parser.add_argument(
        "--betas",
        required=True,
        help="TSV file of normalized beta values",
    )
    parser.add_argument(
        "--pvals",
        required=True,
        help="TSV file of detection p-values",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help="P-value threshold (default: 0.01)",
    )
    parser.add_argument(
        "--sample-name",
        default=None,
        help="Optional explicit sample name for single-sample files",
    )
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    print("=== P-value Correction ===")
    print(f"Beta values file: {args.betas}")
    print(f"Detection p-values file: {args.pvals}")
    if args.sample_name:
        print(f"Sample name: {args.sample_name}")
    print(f"P-value threshold: {args.threshold}\n")

    # Read input files
    beta_df = pd.read_csv(args.betas, sep="\t", index_col=0, header=None)
    pval_df = pd.read_csv(args.pvals, sep="\t", index_col=0, header=None)

    print(
        f"Beta matrix: {beta_df.shape[0]} probes x {beta_df.shape[1]} samples"
    )
    print(
        f"P-value matrix: {pval_df.shape[0]} probes x "
        f"{pval_df.shape[1]} samples\n"
    )

    # Some detection p-value files carry probe suffixes such as `_BC11` or
    # `_TC11` that are not present in the beta matrix. Normalize those rownames
    # to the base probe ID before alignment and collapse any duplicates.
    normalized_pval_index = pval_df.index.str.replace(
        r"_[A-Z0-9]+$", "", regex=True
    )
    if not normalized_pval_index.equals(pval_df.index):
        pval_df = pval_df.copy()
        pval_df.index = normalized_pval_index
        pval_df = pval_df.groupby(level=0).max()

    if args.sample_name and beta_df.shape[1] == 1:
        beta_df = beta_df.copy()
        beta_df.columns = [args.sample_name]

    if args.sample_name and pval_df.shape[1] == 1:
        pval_df = pval_df.copy()
        pval_df.columns = [args.sample_name]

    # Align matrices on shared probes and samples.
    # First remove p-value probes that are not present in the beta matrix.
    beta_probes = beta_df.index
    pval_df = pval_df.loc[pval_df.index.isin(beta_probes)]

    # Preserve the beta-file probe order as the canonical output order.
    common_probes = beta_probes.intersection(pval_df.index)
    common_samples = beta_df.columns.intersection(pval_df.columns)

    if (
        len(common_samples) == 0
        and beta_df.shape[1] == 1
        and pval_df.shape[1] == 1
    ):
        # Per-sample files can legitimately carry different column labels;
        # align the only columns by position.
        print(
            "No shared sample column names found, but both files contain "
            "a single sample. Aligning by position."
        )
        common_samples = beta_df.columns
        pval_df = pval_df.copy()
        pval_df.columns = beta_df.columns

    beta_df = beta_df.loc[common_probes, common_samples]
    pval_df = pval_df.reindex(index=beta_df.index, columns=common_samples)

    if pval_df.isna().all(axis=None):
        raise ValueError(
            "Detection p-values could not be aligned to beta probes. "
            "Check that the beta and p-value files belong to the same sample."
        )

    print(f"Common probes: {len(common_probes)}")
    print(f"Common samples: {len(common_samples)}\n")

    nan_before = beta_df.isna().sum().sum()
    print(f"NaN count before p-value correction: {nan_before}")

    # Replace beta values with NaN where p-value >= threshold
    beta_df = beta_df.where(pval_df < args.threshold, other=np.nan)

    nan_after = beta_df.isna().sum().sum()
    print(f"NaN count after p-value correction: {nan_after}")
    print(f"Values replaced: {nan_after - nan_before}\n")

    # Save corrected beta values
    output_file = f"{args.outdir}/corrected_betas.tsv"
    beta_df.index.name = None
    beta_df.to_csv(output_file, sep="\t", header=False)
    print(f"Saved corrected beta values to {output_file}")

    print("\n=== P-value correction complete ===")


if __name__ == "__main__":
    main()
