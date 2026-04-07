"""
Bulk ATAC-seq Analysis Pipeline
================================
3 conditions × 10 replicates = 30 samples

Steps
-----
1. BAM preprocessing  – deduplication, MAPQ / flag filtering,
                        chrM removal, blacklist masking
2. QC metrics         – flagstat, fragment-size dist, TSS enrichment,
                        FRiP, library complexity, MultiQC report
3. Peak calling       – MACS3 per sample + consensus peak sets
4. Differential analysis – DiffBind (DESeq2) with volcano / MA plots
                           and ChIPseeker annotation

Usage
-----
    snakemake --use-conda --cores 32 -p
    snakemake --use-conda --cores 32 --profile slurm   # HPC
"""

import pandas as pd
from pathlib import Path

# ─── Configuration ─────────────────────────────────────────────────────────────
configfile: "config/config.yaml"

# ─── Sample metadata ───────────────────────────────────────────────────────────
samples_df = pd.read_csv(config["samples"])
samples_df = samples_df.set_index("sample", drop=False)

SAMPLES    = samples_df["sample"].tolist()
CONDITIONS = sorted(samples_df["condition"].unique().tolist())

CONDITION_SAMPLES = {
    cond: samples_df.loc[samples_df["condition"] == cond, "sample"].tolist()
    for cond in CONDITIONS
}

# ─── Global paths ──────────────────────────────────────────────────────────────
RESULTS = config["results_dir"]
BAM_DIR = config.get("bam_dir", "data/bam")

# ─── Include rule modules ──────────────────────────────────────────────────────
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/peak_calling.smk"
include: "workflow/rules/diffbind.smk"

# ─── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        # 1. Final preprocessed BAMs (indexed)
        expand(
            "{results}/bam/final/{sample}.filtered.bam.bai",
            results=RESULTS, sample=SAMPLES,
        ),
        # 2. QC
        expand(
            "{results}/qc/flagstat/{sample}.flagstat",
            results=RESULTS, sample=SAMPLES,
        ),
        expand(
            "{results}/qc/fragment_size/{sample}_fragment_size.txt",
            results=RESULTS, sample=SAMPLES,
        ),
        expand(
            "{results}/qc/frip/{sample}.frip.txt",
            results=RESULTS, sample=SAMPLES,
        ),
        expand(
            "{results}/qc/library_complexity/{sample}.complexity.txt",
            results=RESULTS, sample=SAMPLES,
        ),
        expand(
            "{results}/qc/tss_enrichment/{sample}_TSS_heatmap.pdf",
            results=RESULTS, sample=SAMPLES,
        ),
        f"{RESULTS}/qc/multiqc/multiqc_report.html",
        # 3. Peaks
        expand(
            "{results}/peaks/{sample}/{sample}_peaks.narrowPeak",
            results=RESULTS, sample=SAMPLES,
        ),
        expand(
            "{results}/peaks/consensus/{condition}_consensus_peaks.bed",
            results=RESULTS, condition=CONDITIONS,
        ),
        f"{RESULTS}/peaks/consensus/all_samples_consensus_peaks.bed",
        # 4. DiffBind
        f"{RESULTS}/diffbind/diffbind_analysis.RData",
        f"{RESULTS}/diffbind/plots/pca_plot.pdf",
        f"{RESULTS}/diffbind/plots/correlation_heatmap.pdf",
        f"{RESULTS}/diffbind/differential_peaks_all_contrasts.csv",
