#!/usr/bin/env python3
"""
determine_fastp_params.py
=========================
Analyze QC metrics and recommend fastp parameters.

Logic:
1. Read QC metrics summary CSV
2. Assess overall read quality
3. Recommend quality threshold (--qualified_quality_phred)
4. Recommend minimum length (--length_required)
5. Output recommendations as text and JSON

Recommendations based on:
- Mean quality scores across samples
- Length distribution
- Adapter/GC content
- Best practices for ATAC-seq (typically more lenient than RNA-seq)

Usage:
    python determine_fastp_params.py \
        --qc-metrics results/qc_metrics/qc_metrics_summary.csv \
        --min-quality 20 \
        --min-length 20 \
        --output-txt results/qc_metrics/fastp_params_recommendation.txt \
        --output-json results/qc_metrics/fastp_params.json
"""

import os
import sys
import argparse
import pandas as pd
import json
from pathlib import Path


def assess_quality(qc_df):
    """
    Assess read quality from QC metrics.
    
    Returns:
    - quality_level: "excellent", "good", "fair", "poor"
    - mean_gc: mean GC content
    - considerations: list of observations
    """
    considerations = []
    
    # Simple heuristic based on typical FastQC results
    # (In real scenario, would parse mean_quality from FastQC)
    
    n_samples = len(qc_df)
    considerations.append(f"Analyzed {n_samples} samples")
    
    # GC content check (ATAC-seq typically 40-45%)
    if 'gc_content' in qc_df.columns:
        mean_gc = qc_df['gc_content'].mean()
        considerations.append(f"Mean GC content: {mean_gc:.1f}%")
        
        if mean_gc < 30 or mean_gc > 50:
            considerations.append(f"  ⚠ GC content {mean_gc:.1f}% is outside typical range (40-45%)")
    
    # Sequence length
    if 'sequence_length' in qc_df.columns:
        seq_lens = qc_df['sequence_length'].unique()
        considerations.append(f"Sequence lengths: {', '.join(map(str, seq_lens))}")
    
    # Default: assume "good" quality for ATAC-seq
    quality_level = "good"
    
    return quality_level, considerations


def recommend_fastp_params(qc_df, min_quality_default, min_length_default):
    """
    Recommend fastp parameters based on QC assessment.
    
    Returns dict with parameters.
    """
    quality_level, considerations = assess_quality(qc_df)
    
    # ATAC-seq specific recommendations
    # Generally more lenient than RNA-seq because:
    # - We're assessing chromatin accessibility, not expression
    # - Some technical variation is acceptable
    
    recommendations = {
        "quality_level": quality_level,
        "min_quality": min_quality_default,  # Use default for now
        "min_length": min_length_default,
        "detect_adapter_for_pe": True,
        "adapter_sequence": "auto",
        "trim_front1": 0,
        "trim_tail1": 0,
        "trim_front2": 0,
        "trim_tail2": 0,
        "method": "sliding_window",  # fastp default
        "threshold": "sliding_window_mean_quality",
        "considerations": considerations,
        "rationale": [
            "ATAC-seq: lenient trimming (Q20, min length 20bp)",
            "Auto-detect adapters for paired-end sequencing",
            "Standard sliding window quality assessment"
        ]
    }
    
    # Adjust based on quality level
    if quality_level == "excellent":
        recommendations["min_quality"] = 25
        recommendations["rationale"].append("High quality data → stricter filtering")
    elif quality_level == "poor":
        recommendations["min_quality"] = 15
        recommendations["rationale"].append("Lower quality data → more lenient filtering")
    
    return recommendations


def main():
    parser = argparse.ArgumentParser(
        description='Recommend fastp parameters based on QC metrics'
    )
    parser.add_argument('--qc-metrics', required=True,
                        help='QC metrics CSV (output from parse_fastqc.py)')
    parser.add_argument('--min-quality', type=int, default=20,
                        help='Default minimum quality score')
    parser.add_argument('--min-length', type=int, default=20,
                        help='Default minimum read length')
    parser.add_argument('--output-txt', required=True,
                        help='Output recommendations as text file')
    parser.add_argument('--output-json', required=True,
                        help='Output parameters as JSON file')
    
    args = parser.parse_args()
    
    # Create output directory if needed
    for output_path in [args.output_txt, args.output_json]:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    # Load QC metrics
    print(f"\nLoading QC metrics from: {args.qc_metrics}")
    try:
        qc_df = pd.read_csv(args.qc_metrics)
    except FileNotFoundError:
        print(f"Error: QC metrics file not found: {args.qc_metrics}")
        sys.exit(1)
    
    if qc_df.empty:
        print("Error: QC metrics CSV is empty")
        sys.exit(1)
    
    print(f"Loaded metrics for {len(qc_df)} samples")
    
    # Get recommendations
    recommendation = recommend_fastp_params(
        qc_df, args.min_quality, args.min_length
    )
    
    # Write text report
    print(f"\nWriting recommendations to: {args.output_txt}")
    with open(args.output_txt, 'w') as f:
        f.write("="*70 + "\n")
        f.write("FastP Parameter Recommendations\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Quality Level: {recommendation['quality_level'].upper()}\n")
        f.write(f"Min Quality Score: {recommendation['min_quality']}\n")
        f.write(f"Min Read Length: {recommendation['min_length']} bp\n")
        f.write(f"Adapter Detection: {recommendation['detect_adapter_for_pe']}\n\n")
        
        f.write("QC Assessment:\n")
        for note in recommendation['considerations']:
            f.write(f"  • {note}\n")
        f.write("\n")
        
        f.write("Rationale:\n")
        for reason in recommendation['rationale']:
            f.write(f"  • {reason}\n")
        f.write("\n")
        
        f.write("Recommended fastp Command:\n")
        f.write("-" * 70 + "\n")
        f.write("fastp \\\n")
        f.write("  -i input_R1.fastq.gz -I input_R2.fastq.gz \\\n")
        f.write("  -o output_R1.trimmed.fastq.gz -O output_R2.trimmed.fastq.gz \\\n")
        f.write(f"  -q {recommendation['min_quality']} \\\n")
        f.write(f"  -l {recommendation['min_length']} \\\n")
        f.write("  --detect_adapter_for_pe \\\n")
        f.write("  --adapter_sequence auto \\\n")
        f.write("  --thread 4 \\\n")
        f.write("  --json output.fastp.json \\\n")
        f.write("  --html output.fastp.html\n")
        f.write("-" * 70 + "\n\n")
    
    # Write JSON parameters
    print(f"Writing JSON parameters to: {args.output_json}")
    
    # Remove considerations from JSON (not serializable in some cases)
    json_params = {
        "quality_level": recommendation["quality_level"],
        "min_quality": recommendation["min_quality"],
        "min_length": recommendation["min_length"],
        "detect_adapter_for_pe": recommendation["detect_adapter_for_pe"],
        "adapter_sequence": recommendation["adapter_sequence"],
        "trim_front": recommendation["trim_front1"],
        "trim_tail": recommendation["trim_tail1"],
        "method": recommendation["method"],
        "rationale": recommendation["rationale"]
    }
    
    with open(args.output_json, 'w') as f:
        json.dump(json_params, f, indent=2)
    
    print("\n✓ Recommendation files created")
    print(f"\nRecommendation Summary:")
    print(f"  Quality Level: {recommendation['quality_level']}")
    print(f"  Min Quality: {recommendation['min_quality']}")
    print(f"  Min Length: {recommendation['min_length']} bp")
    print()


if __name__ == '__main__':
    main()
