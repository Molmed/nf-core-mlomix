#!/usr/bin/env python3

"""
P-value correction for methylation beta values.
Replaces beta values with NaN where detection p-value >= threshold.
"""

import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="P-value correction for beta values")
    parser.add_argument("--betas", required=True, help="TSV file of normalized beta values")
    parser.add_argument("--pvals", required=True, help="TSV file of detection p-values")
    parser.add_argument("--threshold", type=float, default=0.01, help="P-value threshold (default: 0.01)")
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    print("=== P-value Correction ===")
    print(f"Beta values file: {args.betas}")
    print(f"Detection p-values file: {args.pvals}")
    print(f"P-value threshold: {args.threshold}\n")

    # Read input files
    beta_df = pd.read_csv(args.betas, sep="\t", index_col=0)
    pval_df = pd.read_csv(args.pvals, sep="\t", index_col=0)

    print(f"Beta matrix: {beta_df.shape[0]} probes x {beta_df.shape[1]} samples")
    print(f"P-value matrix: {pval_df.shape[0]} probes x {pval_df.shape[1]} samples\n")

    # Align matrices on shared probes and samples
    common_probes = beta_df.index.intersection(pval_df.index)
    common_samples = beta_df.columns.intersection(pval_df.columns)
    beta_df = beta_df.loc[common_probes, common_samples]
    pval_df = pval_df.loc[common_probes, common_samples]

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
    beta_df.to_csv(output_file, sep="\t", index_label="probe_id")
    print(f"Saved corrected beta values to {output_file}")

    print("\n=== P-value correction complete ===")


if __name__ == "__main__":
    main()
