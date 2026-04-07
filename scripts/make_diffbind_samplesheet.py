#!/usr/bin/env python3
"""
make_diffbind_samplesheet.py
============================
Generate a DiffBind-compatible CSV samplesheet from the pipeline sample table.

DiffBind samplesheet columns
------------------------------
  SampleID   – unique sample identifier
  Condition  – experimental condition
  Replicate  – replicate number (integer)
  bamReads   – path to filtered BAM
  Peaks      – path to MACS3-filtered narrowPeak file
  PeakCaller – always "macs" for MACS3 output

Usage
-----
  python make_diffbind_samplesheet.py \
      --samples   config/samples.csv \
      --bam-dir   results/bam/final \
      --peaks-dir results/peaks \
      --output    results/diffbind/samplesheet.csv
"""

import argparse
import os
import sys
import pandas as pd


def build_samplesheet(
    samples_csv: str,
    bam_dir: str,
    peaks_dir: str,
    output_csv: str,
) -> None:
    samples_df = pd.read_csv(samples_csv)

    required_cols = {"sample", "condition", "replicate"}
    missing = required_cols - set(samples_df.columns)
    if missing:
        sys.exit(f"ERROR: samples CSV is missing columns: {missing}")

    records = []
    warnings = []

    for _, row in samples_df.iterrows():
        sample    = str(row["sample"])
        condition = str(row["condition"])
        replicate = int(row["replicate"])

        bam_path   = os.path.join(bam_dir,   f"{sample}.filtered.bam")
        peaks_path = os.path.join(peaks_dir, sample, f"{sample}_peaks.filtered.narrowPeak")

        # Warn about missing files but do not abort — files may not exist yet
        # when this script is called by Snakemake to generate the samplesheet.
        for path, label in [(bam_path, "BAM"), (peaks_path, "peaks")]:
            if not os.path.exists(path):
                warnings.append(f"  {label} not found (will be created by pipeline): {path}")

        records.append({
            "SampleID":   sample,
            "Condition":  condition,
            "Replicate":  replicate,
            "bamReads":   bam_path,
            "Peaks":      peaks_path,
            "PeakCaller": "macs",
        })

    out_df = pd.DataFrame(records, columns=[
        "SampleID", "Condition", "Replicate",
        "bamReads", "Peaks", "PeakCaller",
    ])

    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    out_df.to_csv(output_csv, index=False)

    if warnings:
        print("Warnings (non-fatal — paths will exist after pipeline runs):")
        for w in warnings:
            print(w)

    print(f"\nSamplesheet written: {output_csv}")
    print(f"  {len(records)} samples across "
          f"{out_df['Condition'].nunique()} conditions\n")
    print(out_df.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate DiffBind samplesheet from pipeline sample table."
    )
    parser.add_argument("--samples",   required=True, help="Pipeline samples CSV")
    parser.add_argument("--bam-dir",   required=True, help="Directory of filtered BAMs")
    parser.add_argument("--peaks-dir", required=True, help="Root directory of peak files")
    parser.add_argument("--output",    required=True, help="Output DiffBind CSV path")
    args = parser.parse_args()

    build_samplesheet(
        samples_csv = args.samples,
        bam_dir     = args.bam_dir,
        peaks_dir   = args.peaks_dir,
        output_csv  = args.output,
    )


if __name__ == "__main__":
    main()
