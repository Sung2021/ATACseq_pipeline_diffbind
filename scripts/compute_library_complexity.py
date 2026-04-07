#!/usr/bin/env python3
"""
compute_library_complexity.py
=============================
Compute ENCODE ATAC-seq library complexity metrics from a filtered BAM file.

Metrics
-------
  NRF  (Non-Redundant Fraction) = distinct genomic positions / total read pairs
  PBC1 (PCR Bottleneck Coefficient 1) = positions with exactly 1 read / distinct positions
  PBC2 (PCR Bottleneck Coefficient 2) = positions with exactly 1 read / positions with exactly 2 reads

ENCODE quality thresholds
--------------------------
  NRF  >= 0.9   (Concerning < 0.7, Acceptable 0.7–0.9, Ideal >= 0.9)
  PBC1 >= 0.9   (Concerning < 0.7, Acceptable 0.7–0.9, Ideal >= 0.9)
  PBC2 >= 3.0   (Concerning < 1.0, Acceptable 1.0–3.0, Ideal >= 3.0)

Usage
-----
  python compute_library_complexity.py <bam> <sample_name> <output_txt>
"""

import sys
import pysam
from collections import Counter


def qc_flag(value: float, thresholds: tuple) -> str:
    """Return ENCODE quality flag string."""
    lo, hi = thresholds
    if value < lo:
        return "CONCERNING"
    elif value < hi:
        return "ACCEPTABLE"
    return "IDEAL"


def compute_library_complexity(bam_path: str, sample_name: str, output_path: str) -> None:
    position_counts: Counter = Counter()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            # Skip secondary, supplementary, unmapped, and failed QC
            if (read.is_secondary
                    or read.is_supplementary
                    or read.is_unmapped
                    or read.is_qcfail):
                continue
            # For paired-end: count only read1 to avoid double-counting
            if read.is_paired and not read.is_read1:
                continue
            # Skip mitochondrial reads
            ref_name = read.reference_name or ""
            if ref_name in ("chrM", "MT", "chrMT", "M"):
                continue

            strand = "-" if read.is_reverse else "+"
            key = (ref_name, read.reference_start, strand)
            position_counts[key] += 1

    total_reads       = sum(position_counts.values())
    distinct_positions = len(position_counts)
    one_read_pos      = sum(1 for v in position_counts.values() if v == 1)
    two_read_pos      = sum(1 for v in position_counts.values() if v == 2)

    nrf  = distinct_positions / total_reads   if total_reads        > 0 else 0.0
    pbc1 = one_read_pos       / distinct_positions if distinct_positions > 0 else 0.0
    pbc2 = one_read_pos       / two_read_pos  if two_read_pos       > 0 else float("nan")

    nrf_flag  = qc_flag(nrf,  (0.7, 0.9))
    pbc1_flag = qc_flag(pbc1, (0.7, 0.9))
    pbc2_flag = "NA" if two_read_pos == 0 else qc_flag(pbc2, (1.0, 3.0))

    with open(output_path, "w") as fout:
        header = (
            "sample\tTotalReadPairs\tDistinctReadPairs\tOneReadPair\t"
            "TwoReadPairs\tNRF\tPBC1\tPBC2\t"
            "NRF_flag\tPBC1_flag\tPBC2_flag\n"
        )
        fout.write(header)
        pbc2_str = f"{pbc2:.4f}" if two_read_pos > 0 else "NA"
        fout.write(
            f"{sample_name}\t{total_reads}\t{distinct_positions}\t"
            f"{one_read_pos}\t{two_read_pos}\t"
            f"{nrf:.4f}\t{pbc1:.4f}\t{pbc2_str}\t"
            f"{nrf_flag}\t{pbc1_flag}\t{pbc2_flag}\n"
        )

    # Print summary to stdout
    print(f"Sample         : {sample_name}")
    print(f"Total reads    : {total_reads:,}")
    print(f"Distinct pos   : {distinct_positions:,}")
    print(f"NRF            : {nrf:.4f}  [{nrf_flag}]")
    print(f"PBC1           : {pbc1:.4f}  [{pbc1_flag}]")
    print(f"PBC2           : {pbc2_str}  [{pbc2_flag}]")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <bam> <sample_name> <output_txt>",
              file=sys.stderr)
        sys.exit(1)

    compute_library_complexity(sys.argv[1], sys.argv[2], sys.argv[3])
