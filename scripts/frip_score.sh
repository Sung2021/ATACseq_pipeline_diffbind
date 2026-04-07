#!/usr/bin/env bash
# =============================================================================
# frip_score.sh  –  Fraction of Reads in Peaks (FRiP)
#
# FRiP = reads overlapping peaks / total mapped reads
# ENCODE minimum threshold for ATAC-seq: FRiP >= 0.2
#
# Usage:
#   bash frip_score.sh <bam> <peaks.narrowPeak> <sample_name> <output.txt>
#
# Dependencies: samtools, bedtools (both in preprocessing conda env)
# =============================================================================
set -euo pipefail

BAM="${1}"
PEAKS="${2}"
SAMPLE="${3}"
OUTPUT="${4}"

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <bam> <peaks.narrowPeak> <sample_name> <output.txt>" >&2
    exit 1
fi

echo "  FRiP calculation for sample: ${SAMPLE}"

# Total mapped reads (exclude unmapped: -F 4)
total_reads=$(samtools view -c -F 4 "${BAM}")

if [[ "${total_reads}" -eq 0 ]]; then
    echo "WARNING: No mapped reads found in ${BAM}" >&2
    printf "sample\ttotal_reads\treads_in_peaks\tFRiP\tFRiP_flag\n" > "${OUTPUT}"
    printf "%s\t0\t0\tNA\tFAIL\n" "${SAMPLE}" >> "${OUTPUT}"
    exit 0
fi

# Reads overlapping any peak region (-u: count each read once even if
# it overlaps multiple peaks)
reads_in_peaks=$(bedtools intersect \
    -abam "${BAM}" \
    -b    "${PEAKS}" \
    -u \
    | samtools view -c)

# Compute FRiP
frip=$(awk -v rip="${reads_in_peaks}" -v tot="${total_reads}" \
       'BEGIN { printf "%.6f", rip/tot }')

# Quality flag (ENCODE: >= 0.2 acceptable, >= 0.3 ideal for ATAC-seq)
frip_flag=$(awk -v f="${frip}" \
    'BEGIN {
        if (f < 0.1) print "FAIL";
        else if (f < 0.2) print "CONCERNING";
        else if (f < 0.3) print "ACCEPTABLE";
        else print "IDEAL"
    }')

printf "sample\ttotal_reads\treads_in_peaks\tFRiP\tFRiP_flag\n"    > "${OUTPUT}"
printf "%s\t%s\t%s\t%s\t%s\n" \
    "${SAMPLE}" "${total_reads}" "${reads_in_peaks}" "${frip}" "${frip_flag}" \
    >> "${OUTPUT}"

echo "    Total reads    : ${total_reads}"
echo "    Reads in peaks : ${reads_in_peaks}"
echo "    FRiP           : ${frip}  [${frip_flag}]"
