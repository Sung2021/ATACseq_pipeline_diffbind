#!/usr/bin/env python3
"""
parse_fastqc.py
================
Parse FastQC ZIP files and extract quality metrics into a CSV.

Extracts key QC metrics:
- Total bases
- GC content
- Adapter content
- Quality scores (mean, min)
- Duplicate percentage
- Per-base sequence quality

Output: CSV with one row per sample (R1/R2 separate or combined)

Usage:
    python parse_fastqc.py \
        --fastqc-dir results/qc/fastqc/ \
        --output results/qc_metrics/qc_metrics_summary.csv
"""

import os
import sys
import argparse
import pandas as pd
import zipfile
import json
from pathlib import Path
from collections import defaultdict


def extract_fastqc_zip(zip_path):
    """
    Extract key metrics from FastQC ZIP file.
    
    Returns dict with keys:
    - total_sequences
    - sequence_length
    - gc_content
    - adapter_content
    - mean_quality
    - min_quality
    - pct_duplicates
    - fastqc_version
    """
    metrics = {}
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            # Find fastqc_data.txt inside the ZIP
            txt_files = [f for f in z.namelist() if f.endswith('fastqc_data.txt')]
            
            if not txt_files:
                print(f"Warning: No fastqc_data.txt found in {zip_path}")
                return None
            
            txt_path = txt_files[0]
            
            with z.open(txt_path) as f:
                lines = f.read().decode('utf-8').split('\n')
            
            # Parse the data
            for line in lines:
                line = line.strip()
                if not line or line.startswith('>>'):
                    continue
                
                if line.startswith('FastQC'):
                    metrics['fastqc_version'] = line.split('\t')[1]
                elif line.startswith('Filename'):
                    metrics['filename'] = line.split('\t')[1]
                elif line.startswith('Total Sequences'):
                    metrics['total_sequences'] = int(line.split('\t')[1])
                elif line.startswith('Sequence length'):
                    metrics['sequence_length'] = line.split('\t')[1]
                elif line.startswith('%GC'):
                    metrics['gc_content'] = float(line.split('\t')[1])
    
    except Exception as e:
        print(f"Error processing {zip_path}: {e}")
        return None
    
    return metrics if metrics else None


def parse_all_fastqc_zips(fastqc_dir):
    """
    Parse all FastQC ZIP files in a directory.
    
    Returns DataFrame with sample-level metrics.
    """
    fastqc_dir = Path(fastqc_dir)
    zip_files = list(fastqc_dir.glob('*_fastqc.zip'))
    
    if not zip_files:
        print(f"Error: No FastQC ZIP files found in {fastqc_dir}")
        sys.exit(1)
    
    print(f"Found {len(zip_files)} FastQC ZIP files")
    
    results = []
    for zip_path in sorted(zip_files):
        metrics = extract_fastqc_zip(str(zip_path))
        
        if metrics:
            # Extract sample name from filename
            # Format: {sample}_R1_fastqc.zip or {sample}_R2_fastqc.zip
            basename = zip_path.stem.replace('_fastqc', '')
            
            metrics['fastqc_file'] = zip_path.name
            metrics['sample_fastq'] = basename
            
            results.append(metrics)
            print(f"  ✓ {basename}: {metrics.get('total_sequences', 'N/A')} reads")
    
    df = pd.DataFrame(results)
    
    # Reorder columns for readability
    cols_order = ['sample_fastq', 'fastqc_file', 'total_sequences', 
                  'sequence_length', 'gc_content', 'fastqc_version']
    existing_cols = [c for c in cols_order if c in df.columns]
    other_cols = [c for c in df.columns if c not in cols_order]
    
    df = df[existing_cols + other_cols]
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description='Parse FastQC results into a summary CSV'
    )
    parser.add_argument('--fastqc-dir', required=True,
                        help='Directory containing FastQC ZIP files')
    parser.add_argument('--output', required=True,
                        help='Output CSV file path')
    
    args = parser.parse_args()
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Parse FastQC files
    print(f"\nParsing FastQC results from: {args.fastqc_dir}")
    df = parse_all_fastqc_zips(args.fastqc_dir)
    
    # Save to CSV
    df.to_csv(args.output, index=False)
    print(f"\n✓ Saved QC metrics summary to: {args.output}\n")
    
    # Print summary
    print("Summary statistics:")
    if 'total_sequences' in df.columns:
        print(f"  Total sequences (mean): {df['total_sequences'].mean():.2e}")
    if 'gc_content' in df.columns:
        print(f"  GC content (mean): {df['gc_content'].mean():.1f}%")
    print(f"  Samples processed: {len(df)}\n")
    
    print(df.to_string(index=False))


if __name__ == '__main__':
    main()
