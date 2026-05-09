#!/usr/bin/env python3

import os
import sys
import pandas as pd
import argparse

# Fix matplotlib cache directory
os.environ['MPLCONFIGDIR'] = os.path.join(os.getcwd(), '.matplotlib')
os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_class_distribution(y_series, title, outname):
    colormap = plt.colormaps.get_cmap('tab20')
    counts = y_series.value_counts().sort_values(ascending=False)
    plt.figure(figsize=(10, 5))
    counts.plot(kind='bar', color=colormap.colors)

    for i, count in enumerate(counts):
        plt.text(i, count + 0.5, str(count), ha='center', va='bottom', rotation=90)

    plt.title(title)
    plt.xlabel("Class")
    plt.ylabel("Count")
    plt.xticks(rotation=90)
    plt.savefig(outname, bbox_inches='tight')

def plot_empty(message, title, outname):
    plt.figure(figsize=(10, 5))
    plt.text(0.5, 0.5, message, ha='center', va='center', fontsize=14, transform=plt.gca().transAxes)
    plt.title(title)
    plt.xlabel("Class")
    plt.ylabel("Count")
    plt.savefig(outname, bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser(description='Filter samples by class and generate reports')
    parser.add_argument('samplesheet', help='Path to samplesheet file')
    parser.add_argument('--gex-min', type=int, default=3, help='Minimum samples per class for GEX')
    parser.add_argument('--dnam-min', type=int, default=3, help='Minimum samples per class for DNAM')
    parser.add_argument('--gex-class-file', type=str, default=None, help='Optional file with GEX class names to retain (one per line)')
    parser.add_argument('--dnam-class-file', type=str, default=None, help='Optional file with DNAM class names to retain (one per line)')
    args = parser.parse_args()

    # Auto-detect delimiter
    samplesheet_path = args.samplesheet
    if samplesheet_path.endswith('.tsv'):
        sep = '\t'
    else:
        sep = ','

    df = pd.read_csv(samplesheet_path, sep=sep)

    # Determine modality by presence of GEX/DNAM-specific columns
    has_gex = df['gex_feature_counts_file'].notna() if 'gex_feature_counts_file' in df.columns else pd.Series([False] * len(df))
    has_dnam = (df['dnam_beta_matrix_file'].notna() | df['dnam_pvals_file'].notna() | df['idats_basename'].notna()) if any(c in df.columns for c in ['dnam_beta_matrix_file','dnam_pvals_file','idats_basename']) else pd.Series([False] * len(df))

    gex_min = args.gex_min
    dnam_min = args.dnam_min
    gex_class_file = args.gex_class_file
    dnam_class_file = args.dnam_class_file

    def read_class_file(class_file_path):
        """Read class names from file (one per line)."""
        if class_file_path is None:
            return None
        classes = []
        with open(class_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:  # Skip empty lines
                    classes.append(line)
        return set(classes)

    def handle_modality_all_and_filtered(mask, modality_label, tsv_all, svg_all, tsv_filtered, svg_filtered, min_samples, class_file=None):
        subdf = df[mask]
        if 'class' not in subdf.columns:
            print(f"WARNING: Samplesheet does not contain a 'class' column for {modality_label}. Generating empty reports.")
            plot_empty(f"No 'class' column found for {modality_label}", f"Class Distribution ({modality_label})", svg_all)
            plot_empty(f"No 'class' column found for {modality_label}", f"Class Distribution ({modality_label})", svg_filtered)
            pd.DataFrame(columns=["class", "Count"]).to_csv(tsv_all, sep='\t', index=False)
            pd.DataFrame(columns=["class", "Count"]).to_csv(tsv_filtered, sep='\t', index=False)
            return subdf, pd.DataFrame(columns=subdf.columns)

        classes = subdf['class'].dropna()
        if classes.empty:
            print(f"WARNING: 'class' column present for {modality_label} but contains no values. Generating empty reports.")
            plot_empty(f"'class' column has no values for {modality_label}", f"Class Distribution ({modality_label})", svg_all)
            plot_empty(f"'class' column has no values for {modality_label}", f"Class Distribution ({modality_label})", svg_filtered)
            pd.DataFrame(columns=["class", "Count"]).to_csv(tsv_all, sep='\t', index=False)
            pd.DataFrame(columns=["class", "Count"]).to_csv(tsv_filtered, sep='\t', index=False)
            return subdf, pd.DataFrame(columns=subdf.columns)

        # All report
        plot_class_distribution(classes, f"Class Distribution ({modality_label}) - All", svg_all)
        vc_all = classes.value_counts().sort_values(ascending=False)
        vc_all.to_csv(tsv_all, sep='\t', header=["Count"], index_label='class')

        # Determine classes to keep: either from class_file or from min_samples threshold
        if class_file is not None:
            allowed_classes = class_file
            kept_classes = [c for c in vc_all.index if c in allowed_classes]
            if len(kept_classes) == 0:
                print(f"WARNING: No {modality_label} classes from the class file were found in the samplesheet. Filtered report will be empty.")
                plot_empty(f"No {modality_label} classes from file found", f"Class Distribution ({modality_label}) - Filtered", svg_filtered)
                pd.DataFrame(columns=["class", "Count"]).to_csv(tsv_filtered, sep='\t', index=False)
                return subdf, pd.DataFrame(columns=subdf.columns)
        else:
            kept_classes = vc_all[vc_all >= min_samples].index.tolist()
            if len(kept_classes) == 0:
                print(f"WARNING: No {modality_label} classes meet the minimum sample threshold of {min_samples}. Filtered report will be empty.")
                plot_empty(f"No {modality_label} classes meet the threshold", f"Class Distribution ({modality_label}) - Filtered", svg_filtered)
                pd.DataFrame(columns=["class", "Count"]).to_csv(tsv_filtered, sep='\t', index=False)
                return subdf, pd.DataFrame(columns=subdf.columns)

        filtered_subdf = subdf[subdf['class'].isin(kept_classes)]
        plot_class_distribution(filtered_subdf['class'], f"Class Distribution ({modality_label}) - Filtered", svg_filtered)
        vc_f = filtered_subdf['class'].value_counts().sort_values(ascending=False)
        vc_f.to_csv(tsv_filtered, sep='\t', header=["Count"], index_label='class')
        return subdf, filtered_subdf

    # Handle GEX
    gex_class_file_set = read_class_file(gex_class_file)
    if has_gex.any():
        gex_all_df, gex_filtered_df = handle_modality_all_and_filtered(has_gex, 'GEX', 'classes.gex.all.tsv', 'class_distribution.gex.all.svg', 'classes.gex.filtered.tsv', 'class_distribution.gex.filtered.svg', gex_min, class_file=gex_class_file_set)
    else:
        plot_empty('No GEX samples found', 'Class Distribution (GEX) - All', 'class_distribution.gex.all.svg')
        plot_empty('No GEX samples found', 'Class Distribution (GEX) - Filtered', 'class_distribution.gex.filtered.svg')
        pd.DataFrame(columns=["class", "Count"]).to_csv('classes.gex.all.tsv', sep='\t', index=False)
        pd.DataFrame(columns=["class", "Count"]).to_csv('classes.gex.filtered.tsv', sep='\t', index=False)
        gex_all_df = pd.DataFrame(columns=df.columns)
        gex_filtered_df = pd.DataFrame(columns=df.columns)

    # Handle DNAM
    dnam_class_file_set = read_class_file(dnam_class_file)
    if has_dnam.any():
        dnam_all_df, dnam_filtered_df = handle_modality_all_and_filtered(has_dnam, 'DNAM', 'classes.dnam.all.tsv', 'class_distribution.dnam.all.svg', 'classes.dnam.filtered.tsv', 'class_distribution.dnam.filtered.svg', dnam_min, class_file=dnam_class_file_set)
    else:
        plot_empty('No DNAM samples found', 'Class Distribution (DNAM) - All', 'class_distribution.dnam.all.svg')
        plot_empty('No DNAM samples found', 'Class Distribution (DNAM) - Filtered', 'class_distribution.dnam.filtered.svg')
        pd.DataFrame(columns=["class", "Count"]).to_csv('classes.dnam.all.tsv', sep='\t', index=False)
        pd.DataFrame(columns=["class", "Count"]).to_csv('classes.dnam.filtered.tsv', sep='\t', index=False)
        dnam_all_df = pd.DataFrame(columns=df.columns)
        dnam_filtered_df = pd.DataFrame(columns=df.columns)

    # Build combined filtered samplesheet: keep samples that meet the minimum for all modalities they participate in
    keep_mask = pd.Series([True] * len(df))
    # If sample has GEX, require it to be in gex_filtered_df
    if not gex_filtered_df.empty and 'sample' in gex_filtered_df.columns:
        gex_ids = set(gex_filtered_df['sample'].astype(str).tolist())
    else:
        gex_ids = set()
    if not dnam_filtered_df.empty and 'sample' in dnam_filtered_df.columns:
        dnam_ids = set(dnam_filtered_df['sample'].astype(str).tolist())
    else:
        dnam_ids = set()

    for idx, row in df.iterrows():
        sid = str(row['sample']) if 'sample' in df.columns else None
        if sid is None:
            keep_mask.iloc[idx] = False
            continue
        # if sample has GEX data, it must be in gex_ids
        if ('gex_feature_counts_file' in df.columns) and pd.notna(row.get('gex_feature_counts_file')):
            if sid not in gex_ids:
                keep_mask.iloc[idx] = False
                continue
        # if sample has DNAM data, it must be in dnam_ids
        if any(c in df.columns for c in ['dnam_beta_matrix_file','dnam_pvals_file','idats_basename']) and (pd.notna(row.get('dnam_beta_matrix_file')) or pd.notna(row.get('dnam_pvals_file')) or pd.notna(row.get('idats_basename'))):
            if sid not in dnam_ids:
                keep_mask.iloc[idx] = False
                continue

    filtered_df = df[keep_mask]

    # Sort by sample column to keep output consistent
    if 'sample' in filtered_df.columns:
        filtered_df = filtered_df.sort_values(by='sample')

    # Save filtered samplesheet as TSV
    filtered_df.to_csv('samplesheet.filtered.csv', sep=',', index=False)

    # Create versions file
    import numpy
    with open("versions.yml", "w") as f:
        f.write(f'"{os.path.basename(__file__)}":\n')
        f.write(f'    python: "{sys.version.split()[0]}"\n')
        f.write(f'    pandas: "{pd.__version__}"\n')
        f.write(f'    numpy: "{numpy.__version__}"\n')
        f.write(f'    matplotlib: "{matplotlib.__version__}"\n')

if __name__ == '__main__':
    main()
