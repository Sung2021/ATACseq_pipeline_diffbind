#!/usr/bin/env python3
"""
compare_qc_before_after.py
===========================
Compare QC metrics before and after fastp trimming.

Extracts and compares:
- Read count (how many reads removed)
- Quality scores (mean, distribution)
- GC content
- Adapter content (should decrease after trimming)
- Duplicate rate

Output:
- comparison_summary.csv: Before/After metrics per sample
- comparison_plots.html: Interactive plots for comparison
- improvement_metrics.txt: Summary of improvements

Usage:
    python compare_qc_before_after.py \
        --fastqc-raw results/qc/fastqc/ \
        --fastqc-trimmed results/qc/fastqc/ \
        --fastp-reports results/trimmed/ \
        --output-csv results/qc_metrics/qc_comparison_summary.csv \
        --output-html results/qc_metrics/qc_comparison_report.html \
        --output-txt results/qc_metrics/improvement_summary.txt
"""

import os
import sys
import argparse
import pandas as pd
import zipfile
import json
from pathlib import Path
from collections import defaultdict
import re


def extract_qc_metrics_from_zip(zip_path):
    """Extract and parse FastQC ZIP file for key metrics."""
    metrics = {}
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            # Find fastqc_data.txt
            txt_files = [f for f in z.namelist() if f.endswith('fastqc_data.txt')]
            if not txt_files:
                return None
            
            with z.open(txt_files[0]) as f:
                content = f.read().decode('utf-8')
            
            # Parse basic metrics
            for line in content.split('\n'):
                line = line.strip()
                if not line or line.startswith('>>'):
                    continue
                
                if line.startswith('Total Sequences'):
                    metrics['total_sequences'] = int(line.split('\t')[1])
                elif line.startswith('Sequence length'):
                    metrics['sequence_length'] = line.split('\t')[1]
                elif line.startswith('%GC'):
                    metrics['gc_content'] = float(line.split('\t')[1])
            
            # Parse per-sequence quality (more detailed)
            if 'Per sequence quality scores' in content:
                # Extract mean quality info from the module
                lines = content.split('\n')
                in_quality_section = False
                qualities = []
                
                for i, line in enumerate(lines):
                    if 'Per sequence quality scores' in line:
                        in_quality_section = True
                        continue
                    if in_quality_section and line.startswith('>>'):
                        break
                    if in_quality_section and '\t' in line and not line.startswith('#'):
                        try:
                            parts = line.split('\t')
                            if len(parts) >= 2:
                                q_score = float(parts[0])
                                count = int(parts[1])
                                # Weight by count to get average
                                qualities.extend([q_score] * min(count, 10))  # Limit to avoid memory issues
                        except:
                            pass
                
                if qualities:
                    metrics['mean_quality'] = sum(qualities) / len(qualities)
                    metrics['min_quality'] = min(qualities)
                    metrics['max_quality'] = max(qualities)
    
    except Exception as e:
        print(f"Warning: Error parsing {zip_path}: {e}")
        return None
    
    return metrics if metrics else None


def extract_fastp_metrics(json_path):
    """Extract metrics from fastp JSON report."""
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # Extract key metrics
        metrics = {
            'input_total_reads': data.get('summary', {}).get('before_filtering', {}).get('total_reads', 0),
            'output_total_reads': data.get('summary', {}).get('after_filtering', {}).get('total_reads', 0),
            'reads_removed': data.get('summary', {}).get('before_filtering', {}).get('total_reads', 0) - 
                           data.get('summary', {}).get('after_filtering', {}).get('total_reads', 0),
            'adapter_trimmed': data.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0),
            'quality_trimmed': len([x for x in data.get('filtering_result', []) if x.get('type') == 'Q20']),
        }
        return metrics
    except Exception as e:
        print(f"Warning: Error parsing {json_path}: {e}")
        return None


def parse_all_comparisons(fastqc_raw_dir, fastqc_trimmed_dir, fastp_dir):
    """
    Parse all QC files and create comparison table.
    
    Returns DataFrame with before/after metrics.
    """
    fastqc_raw_dir = Path(fastqc_raw_dir)
    fastqc_trimmed_dir = Path(fastqc_trimmed_dir)
    fastp_dir = Path(fastp_dir)
    
    # Get raw FastQC files (not trimmed)
    raw_zips = sorted([z for z in fastqc_raw_dir.glob('*_fastqc.zip') 
                       if 'trimmed' not in z.name])
    
    if not raw_zips:
        print(f"No raw FastQC files found in {fastqc_raw_dir}")
        return None
    
    print(f"Found {len(raw_zips)} raw FastQC files")
    
    results = []
    
    for raw_zip in raw_zips:
        # Extract sample name
        basename = raw_zip.stem.replace('_fastqc', '')
        sample_name = basename.replace('_R1', '').replace('_R2', '')
        read_type = 'R1' if '_R1' in basename else 'R2' if '_R2' in basename else 'unknown'
        
        # Parse raw QC
        raw_metrics = extract_qc_metrics_from_zip(str(raw_zip))
        if not raw_metrics:
            continue
        
        # Find trimmed counterpart
        trimmed_zip = fastqc_trimmed_dir / f"{basename}.trimmed_fastqc.zip"
        if not trimmed_zip.exists():
            print(f"Warning: No trimmed FastQC found for {basename}")
            continue
        
        # Parse trimmed QC
        trimmed_metrics = extract_qc_metrics_from_zip(str(trimmed_zip))
        if not trimmed_metrics:
            continue
        
        # Find fastp JSON report
        fastp_json = fastp_dir / f"{sample_name}.fastp.json"
        fastp_metrics = {}
        if fastp_json.exists():
            fastp_metrics = extract_fastp_metrics(str(fastp_json))
        
        # Calculate improvements
        improvement = {
            'sample': sample_name,
            'read_type': read_type,
            'raw_sequences': raw_metrics.get('total_sequences', 0),
            'trimmed_sequences': trimmed_metrics.get('total_sequences', 0),
            'reads_removed': raw_metrics.get('total_sequences', 0) - trimmed_metrics.get('total_sequences', 0),
            'removal_rate_pct': 100 * (raw_metrics.get('total_sequences', 0) - trimmed_metrics.get('total_sequences', 0)) / max(raw_metrics.get('total_sequences', 1), 1),
            'raw_gc': raw_metrics.get('gc_content', None),
            'trimmed_gc': trimmed_metrics.get('gc_content', None),
            'raw_quality': raw_metrics.get('mean_quality', None),
            'trimmed_quality': trimmed_metrics.get('mean_quality', None),
            'quality_improvement': (trimmed_metrics.get('mean_quality', 0) - raw_metrics.get('mean_quality', 0)) if raw_metrics.get('mean_quality') else None,
        }
        
        results.append(improvement)
        
        print(f"  ✓ {basename}: {improvement['reads_removed']:,} reads removed "
              f"({improvement['removal_rate_pct']:.1f}%)")
    
    df = pd.DataFrame(results)
    return df


def generate_html_report(df, output_html):
    """Generate interactive HTML report with plots."""
    html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>FASTQ QC Comparison Report</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.26.0/plotly.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; border-radius: 8px; }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        .metric { display: inline-block; padding: 15px; margin: 10px; background-color: #ecf0f1; border-radius: 5px; }
        .metric-value { font-size: 24px; font-weight: bold; color: #27ae60; }
        .metric-label { font-size: 12px; color: #7f8c8d; }
        .plot { margin: 20px 0; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { padding: 10px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #3498db; color: white; }
        tr:hover { background-color: #f9f9f9; }
        .warning { background-color: #fff3cd; padding: 10px; border-left: 4px solid #ffc107; margin: 20px 0; }
        .success { background-color: #d4edda; padding: 10px; border-left: 4px solid #28a745; margin: 20px 0; }
    </style>
</head>
<body>
    <div class="container">
        <h1>📊 FASTQ QC Comparison Report (Before / After Trimming)</h1>
        
        <div class="success">
            <strong>✓ Analysis Complete:</strong> Comparison of {n_samples} samples across R1/R2 reads
        </div>
        
        <h2>📈 Summary Metrics</h2>
        <div>
            <div class="metric">
                <div class="metric-value">{total_reads_removed:,}</div>
                <div class="metric-label">Total Reads Removed</div>
            </div>
            <div class="metric">
                <div class="metric-value">{avg_removal_rate:.1f}%</div>
                <div class="metric-label">Average Removal Rate</div>
            </div>
            <div class="metric">
                <div class="metric-value">{avg_quality_improvement:.2f}</div>
                <div class="metric-label">Avg Quality Score Improvement</div>
            </div>
        </div>
        
        <h2>📋 Detailed Comparison Table</h2>
        {table_html}
        
        <h2>📉 Visualization</h2>
        <div id="plot1" class="plot" style="width:100%; height:500px;"></div>
        <div id="plot2" class="plot" style="width:100%; height:500px;"></div>
        <div id="plot3" class="plot" style="width:100%; height:500px;"></div>
        
        <h2>💡 Interpretation</h2>
        <ul>
            <li><strong>Removal Rate:</strong> Percentage of reads filtered out
                <ul>
                    <li>5-15%: Typical for good quality data</li>
                    <li>15-30%: Data with moderate quality issues</li>
                    <li>&gt;30%: Consider checking raw data quality</li>
                </ul>
            </li>
            <li><strong>Quality Improvement:</strong> Change in mean quality score
                <ul>
                    <li>Positive: Quality improved after trimming ✓</li>
                    <li>Negative: Quality decreased (unusual) ⚠️</li>
                </ul>
            </li>
            <li><strong>GC Content:</strong> Should remain relatively stable</li>
        </ul>
        
        <h2>🎯 Recommendations</h2>
        {recommendations}
    </div>
    
    <script>
        // Plot 1: Reads removed per sample
        var samples = {samples_list};
        var removal_pct = {removal_pct_list};
        
        var trace1 = {{
            x: samples,
            y: removal_pct,
            type: 'bar',
            marker: {{color: 'rgba(52, 152, 219, 0.8)'}}
        }};
        
        var layout1 = {{
            title: 'Reads Removed by Trimming (%)',
            xaxis: {{title: 'Sample'}},
            yaxis: {{title: 'Removal Rate (%)'}},
            hovermode: 'closest'
        }};
        
        Plotly.newPlot('plot1', [trace1], layout1, {{responsive: true}});
        
        // Plot 2: Quality improvement
        var quality_improvement = {quality_improvement_list};
        
        var trace2 = {{
            x: samples,
            y: quality_improvement,
            type: 'bar',
            marker: {{
                color: quality_improvement.map(x => x >= 0 ? 'green' : 'red')
            }}
        }};
        
        var layout2 = {{
            title: 'Quality Score Improvement (Before → After)',
            xaxis: {{title: 'Sample'}},
            yaxis: {{title: 'Δ Mean Quality Score'}},
            hovermode: 'closest'
        }};
        
        Plotly.newPlot('plot2', [trace2], layout2, {{responsive: true}});
        
        // Plot 3: Before/After reads comparison
        var before_reads = {before_reads_list};
        var after_reads = {after_reads_list};
        
        var trace3_before = {{
            x: samples,
            y: before_reads,
            name: 'Before Trimming',
            type: 'bar',
            marker: {{color: 'rgba(192, 57, 43, 0.8)'}}
        }};
        
        var trace3_after = {{
            x: samples,
            y: after_reads,
            name: 'After Trimming',
            type: 'bar',
            marker: {{color: 'rgba(39, 174, 96, 0.8)'}}
        }};
        
        var layout3 = {{
            title: 'Read Count Comparison',
            xaxis: {{title: 'Sample'}},
            yaxis: {{title: 'Number of Reads'}},
            barmode: 'group',
            hovermode: 'closest'
        }};
        
        Plotly.newPlot('plot3', [trace3_before, trace3_after], layout3, {{responsive: true}});
    </script>
</body>
</html>
"""
    
    # Prepare data for the report
    n_samples = len(df)
    total_reads_removed = df['reads_removed'].sum()
    avg_removal_rate = df['removal_rate_pct'].mean()
    avg_quality_improvement = df['quality_improvement'].mean() if 'quality_improvement' in df.columns else 0
    
    # Create table
    table_html = df.to_html(index=False, border=0, 
                            float_format=lambda x: f'{x:.2f}' if isinstance(x, float) else str(x))
    
    # Generate recommendations
    recommendations = "<ul>"
    if avg_removal_rate > 30:
        recommendations += "<li>⚠️ High removal rate (&gt;30%). Consider checking raw data quality or adjusting parameters.</li>"
    elif avg_removal_rate < 5:
        recommendations += "<li>ℹ️ Low removal rate (&lt;5%). Data quality is excellent.</li>"
    else:
        recommendations += "<li>✓ Removal rate is within typical range (5-30%).</li>"
    
    if avg_quality_improvement > 0:
        recommendations += "<li>✓ Quality scores improved after trimming (good sign).</li>"
    elif avg_quality_improvement < 0:
        recommendations += "<li>⚠️ Quality scores decreased. This is unusual; check parameters.</li>"
    else:
        recommendations += "<li>ℹ️ Quality scores remained unchanged.</li>"
    
    recommendations += "<li><strong>Next Step:</strong> Review before/after plots. If satisfied, use trimmed reads for next analysis step (e.g., alignment).</li>"
    recommendations += "</ul>"
    
    # Prepare plot data
    samples_list = df['sample'].unique().tolist()
    removal_pct_list = df['removal_rate_pct'].tolist()
    quality_improvement_list = df['quality_improvement'].fillna(0).tolist()
    before_reads_list = df['raw_sequences'].tolist()
    after_reads_list = df['trimmed_sequences'].tolist()
    
    # Format the HTML
    html_final = html_content.format(
        n_samples=n_samples,
        total_reads_removed=int(total_reads_removed),
        avg_removal_rate=avg_removal_rate,
        avg_quality_improvement=avg_quality_improvement,
        table_html=table_html,
        recommendations=recommendations,
        samples_list=samples_list,
        removal_pct_list=removal_pct_list,
        quality_improvement_list=quality_improvement_list,
        before_reads_list=before_reads_list,
        after_reads_list=after_reads_list
    )
    
    with open(output_html, 'w') as f:
        f.write(html_final)


def main():
    parser = argparse.ArgumentParser(
        description='Compare QC metrics before and after fastp trimming'
    )
    parser.add_argument('--fastqc-raw', required=True,
                        help='Directory with raw FastQC ZIP files')
    parser.add_argument('--fastqc-trimmed', required=True,
                        help='Directory with trimmed FastQC ZIP files')
    parser.add_argument('--fastp-reports', required=True,
                        help='Directory with fastp JSON reports')
    parser.add_argument('--output-csv', required=True,
                        help='Output CSV file')
    parser.add_argument('--output-html', required=True,
                        help='Output HTML report')
    parser.add_argument('--output-txt', required=True,
                        help='Output text summary')
    
    args = parser.parse_args()
    
    # Create output directory
    for output_path in [args.output_csv, args.output_html, args.output_txt]:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    print("\nComparing QC metrics before and after trimming...")
    
    # Parse comparisons
    df = parse_all_comparisons(args.fastqc_raw, args.fastqc_trimmed, args.fastp_reports)
    
    if df is None or df.empty:
        print("Error: No comparisons found")
        sys.exit(1)
    
    # Save CSV
    print(f"\nSaving comparison CSV to: {args.output_csv}")
    df.to_csv(args.output_csv, index=False)
    
    # Generate HTML report
    print(f"Generating HTML report: {args.output_html}")
    generate_html_report(df, args.output_html)
    
    # Write text summary
    print(f"Writing text summary to: {args.output_txt}")
    with open(args.output_txt, 'w') as f:
        f.write("="*70 + "\n")
        f.write("QC Comparison Summary: Before vs After Trimming\n")
        f.write("="*70 + "\n\n")
        
        f.write("Overall Statistics:\n")
        f.write(f"  Total Reads Removed: {df['reads_removed'].sum():,}\n")
        f.write(f"  Average Removal Rate: {df['removal_rate_pct'].mean():.1f}%\n")
        if 'quality_improvement' in df.columns:
            f.write(f"  Average Quality Improvement: {df['quality_improvement'].mean():.2f}\n")
        f.write("\n")
        
        f.write("Per-Sample Breakdown:\n")
        f.write("-" * 70 + "\n")
        for _, row in df.iterrows():
            f.write(f"\n{row['sample']} ({row['read_type']}):\n")
            f.write(f"  Reads: {row['raw_sequences']:,} → {row['trimmed_sequences']:,}\n")
            f.write(f"  Removed: {row['reads_removed']:,} ({row['removal_rate_pct']:.1f}%)\n")
            if pd.notna(row.get('quality_improvement')):
                f.write(f"  Quality: {row.get('raw_quality', 'N/A')} → {row.get('trimmed_quality', 'N/A')} "
                       f"(Δ {row['quality_improvement']:.2f})\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("Interpretation:\n")
        f.write("  • Removal Rate <5%: Excellent raw quality\n")
        f.write("  • Removal Rate 5-15%: Good quality (typical)\n")
        f.write("  • Removal Rate 15-30%: Moderate quality issues\n")
        f.write("  • Removal Rate >30%: Consider re-trimming or checking parameters\n")
        f.write("\n")
        f.write("Decision:\n")
        avg_removal = df['removal_rate_pct'].mean()
        if avg_removal > 30:
            f.write("  ⚠️  High removal rate. You MAY want to:\n")
            f.write("      1. Try stricter quality thresholds\n")
            f.write("      2. Increase trimming length\n")
            f.write("      3. OR accept quality and continue (decision is yours)\n")
        elif avg_removal > 15:
            f.write("  ℹ️  Moderate removal rate. Trimming parameters are reasonable.\n")
            f.write("      Can proceed with next step.\n")
        else:
            f.write("  ✓ Low removal rate. Data quality is good. Proceed with next step.\n")
    
    print("\n✓ Comparison complete!")
    print(f"\nOutputs:")
    print(f"  CSV: {args.output_csv}")
    print(f"  HTML: {args.output_html}")
    print(f"  Summary: {args.output_txt}")
    print()


if __name__ == '__main__':
    main()
